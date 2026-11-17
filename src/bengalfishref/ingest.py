"""Rate-limited sequence ingestion via official ncbi_fetch.py wrapper and BOLD v4 Public API."""

from __future__ import annotations

from datetime import datetime, timezone
import importlib.util
import json
from pathlib import Path
import re
import sqlite3
import sys
import time
from typing import Any

import requests

from bengalfishref.audit import evaluate_sequence_hygiene, make_seq_uid
from bengalfishref.taxonomy import normalize_binomial, resolve_species_name_Map

NCBI_FETCH_SCRIPT = Path(
    "/home/sabuj/.gemini/config/plugins/science/skills/ncbi_sequence_fetch/scripts/ncbi_fetch.py"
)

MARKER_ENTREZ_FILTERS: dict[str, str] = {
    "COI": '(COI[Gene] OR COX1[Gene] OR "cytochrome oxidase subunit I"[Title] OR "cytochrome c oxidase subunit I"[Title]) AND 480:1700[SLEN]',
    "12S": '("12S"[Title] OR "12S ribosomal RNA"[Title] OR "s-rRNA"[Gene]) AND 140:1200[SLEN]',
    "16S": '("16S"[Title] OR "16S ribosomal RNA"[Title] OR "l-rRNA"[Gene]) AND 160:1900[SLEN]',
    "Cytb": '(CYTB[Gene] OR "cytochrome b"[Title]) AND 450:1350[SLEN]',
    "Mitogenome": '(mitochondrion[Filter] AND ("complete genome"[Title] OR "mitochondrial DNA, complete"[Title])) AND 14000:22000[SLEN]',
}


def extract_accession_and_species_from_header(header: str) -> tuple[str, str, str]:
    """Extract (accession, reported_binomial, country_hint) from an NCBI FASTA header."""
    hdr = header.lstrip(">").strip()
    parts = hdr.split(None, 1)
    acc = parts[0] if parts else "UNKNOWN"
    desc = parts[1] if len(parts) > 1 else ""
    # Strip leading prefixes like UNVERIFIED: or PREDICTED:
    desc_clean = re.sub(r"^(UNVERIFIED:\s*|PREDICTED:\s*|TPA:\s*)", "", desc, flags=re.IGNORECASE)
    words = desc_clean.split()
    binomial = normalize_binomial(" ".join(words[:2])) if len(words) >= 2 else ""

    country = "Global"
    for c_tag in ("Bangladesh", "Bay of Bengal", "Sundarbans", "Padma", "Meghna", "India", "West Bengal", "Myanmar"):
        if c_tag.lower() in desc.lower():
            country = "Bangladesh" if c_tag in ("Bangladesh", "Sundarbans", "Padma", "Meghna") else c_tag
            break
    return acc, binomial, country


_NCBI_MOD: Any = None


def _get_ncbi_module() -> Any:
    """Load the official science skill `ncbi_fetch.py` module in-process so its polite_http rate limiter persists."""
    global _NCBI_MOD
    if _NCBI_MOD is None:
        if not NCBI_FETCH_SCRIPT.exists():
            raise FileNotFoundError(f"Official wrapper not found at {NCBI_FETCH_SCRIPT}")
        spec = importlib.util.spec_from_file_location("ncbi_fetch_official", str(NCBI_FETCH_SCRIPT))
        if spec is None or spec.loader is None:
            raise RuntimeError(f"Unable to load spec for {NCBI_FETCH_SCRIPT}")
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        _NCBI_MOD = mod
    return _NCBI_MOD


def run_ncbi_fetch_batch(query: str, retmax: int = 120) -> tuple[list[dict[str, Any]], int]:
    """Invoke official `ncbi_fetch.py` esearch + efetch functions with built-in 10 req/s rate limiting and retry."""
    mod = _get_ncbi_module()
    for attempt in range(3):
        try:
            ids, total_count = mod.esearch("nuccore", query, retmax=retmax)
            sequences: list[dict[str, Any]] = []
            if ids:
                batch_size = 40
                for i in range(0, len(ids), batch_size):
                    sub_ids = ids[i : i + batch_size]
                    fasta = mod.efetch("nuccore", ",".join(sub_ids))
                    if fasta:
                        for header, seq in mod.parse_fasta(fasta):
                            sequences.append({"header": header, "sequence": seq, "length": len(seq)})
            return sequences, int(total_count)
        except Exception as exc:
            if attempt == 2:
                print(f"[WARN] NCBI batch failed after 3 attempts: {exc}", file=sys.stderr, flush=True)
                return [], 0
            time.sleep(2.0 * (attempt + 1))
    return [], 0


def insert_ingested_sequence(
    conn: sqlite3.Connection,
    *,
    source_db: str,
    accession: str,
    reported_taxon: str,
    taxon_meta: dict[str, Any],
    marker: str,
    raw_seq: str,
    country_origin: str = "Global",
) -> bool:
    """Evaluate hygiene and upsert a sequence record into SQLite."""
    is_bd = 1 if country_origin == "Bangladesh" else 0
    hygiene = evaluate_sequence_hygiene(
        raw_seq,
        marker=marker,
        genetic_code=int(taxon_meta["genetic_code"]),
        is_voucher_or_geo=(country_origin in ("Bangladesh", "Bay of Bengal", "India")),
    )
    if hygiene["seq_length"] == 0:
        return False

    uid = make_seq_uid(source_db, accession, marker)
    conn.execute(
        """
        INSERT INTO sequences (
            seq_uid, source_db, accession, reported_taxon, resolved_taxon_id,
            marker, seq_length, ambig_pct, gc_pct, country_origin,
            is_bangladesh_voucher, has_stop_codon_numt, homopolymer_max,
            quality_tier, quarantine_reason, sequence
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT(seq_uid) DO UPDATE SET
            quality_tier = excluded.quality_tier,
            quarantine_reason = excluded.quarantine_reason,
            has_stop_codon_numt = excluded.has_stop_codon_numt
        """,
        (
            uid,
            source_db,
            accession,
            reported_taxon,
            int(taxon_meta["taxon_id"]),
            marker,
            hygiene["seq_length"],
            hygiene["ambig_pct"],
            hygiene["gc_pct"],
            country_origin,
            is_bd,
            hygiene["has_stop_codon_numt"],
            hygiene["homopolymer_max"],
            hygiene["quality_tier"],
            hygiene["quarantine_reason"],
            hygiene["clean_seq"],
        ),
    )
    return True


def ingest_ncbi_for_checklist(
    conn: sqlite3.Connection,
    markers: tuple[str, ...] = ("COI", "12S", "16S", "Cytb", "Mitogenome"),
    species_batch_size: int = 8,
    retmax_per_batch: int = 120,
    max_batches: int | None = None,
) -> dict[str, int]:
    """Resumably ingest NCBI nuccore sequences across all checklist species and synonyms."""
    name_map = resolve_species_name_Map(conn)
    taxa_rows = [dict(r) for r in conn.execute("SELECT * FROM taxa ORDER BY valid_name").fetchall()]
    done_batches = {
        r["batch_key"]
        for r in conn.execute("SELECT batch_key FROM ingestion_checkpoints WHERE source_db='NCBI'").fetchall()
    }

    stats = {"batches_executed": 0, "sequences_stored": 0, "batches_skipped_checkpoint": 0}

    for marker in markers:
        m_filter = MARKER_ENTREZ_FILTERS[marker]
        for b_idx in range(0, len(taxa_rows), species_batch_size):
            batch_taxa = taxa_rows[b_idx : b_idx + species_batch_size]
            batch_key = f"NCBI::{marker}::batch_{b_idx // species_batch_size:04d}_n{len(batch_taxa)}"
            if batch_key in done_batches:
                stats["batches_skipped_checkpoint"] += 1
                continue
            if max_batches is not None and stats["batches_executed"] >= max_batches:
                return stats

            org_clauses: list[str] = []
            for tx in batch_taxa:
                org_clauses.append(f'"{tx["valid_name"]}"[Organism]')
                for syn in json.loads(tx["synonyms_json"] or "[]"):
                    ns = normalize_binomial(syn)
                    if ns:
                        org_clauses.append(f'"{ns}"[Organism]')

            entrez_query = f"({' OR '.join(org_clauses)}) AND ({m_filter}) NOT patent[Properties]"
            seqs, total_count = run_ncbi_fetch_batch(entrez_query, retmax=retmax_per_batch)

            # Adaptive Split-on-Saturation: if a hyper-sequenced species (e.g., Danio rerio / Oreochromis niloticus)
            # saturates retmax_per_batch, query each species in this batch individually (retmax=25) so zero
            # co-batched Bangladeshi species suffer false Dark-Taxa dropout!
            if total_count > retmax_per_batch and len(batch_taxa) > 1:
                seqs = []
                per_sp_cap = 10 if marker == "Mitogenome" else 25
                for tx in batch_taxa:
                    sp_clauses = [f'"{tx["valid_name"]}"[Organism]']
                    for syn in json.loads(tx["synonyms_json"] or "[]"):
                        ns = normalize_binomial(syn)
                        if ns:
                            sp_clauses.append(f'"{ns}"[Organism]')
                    sp_q = f"({' OR '.join(sp_clauses)}) AND ({m_filter}) NOT patent[Properties]"
                    sp_seqs, _ = run_ncbi_fetch_batch(sp_q, retmax=per_sp_cap)
                    seqs.extend(sp_seqs)
            n_saved = 0
            with conn:
                for item in seqs:
                    acc, rep_bin, country = extract_accession_and_species_from_header(item.get("header", ""))
                    tx_meta = name_map.get(rep_bin.lower())
                    if tx_meta is None:
                        # Try matching any checklist name appearing inside header string
                        hdr_low = item.get("header", "").lower()
                        for cand_name, cand_meta in name_map.items():
                            if cand_name in hdr_low:
                                tx_meta = cand_meta
                                rep_bin = cand_meta["valid_name"]
                                break
                    if tx_meta is None:
                        continue
                    if insert_ingested_sequence(
                        conn,
                        source_db="NCBI",
                        accession=acc,
                        reported_taxon=rep_bin,
                        taxon_meta=tx_meta,
                        marker=marker,
                        raw_seq=item.get("sequence", ""),
                        country_origin=country,
                    ):
                        n_saved += 1

                conn.execute(
                    """
                    INSERT OR REPLACE INTO ingestion_checkpoints (batch_key, source_db, marker, completed_at, records_fetched)
                    VALUES (?, 'NCBI', ?, ?, ?)
                    """,
                    (batch_key, marker, datetime.now(timezone.utc).isoformat(), n_saved),
                )

            stats["batches_executed"] += 1
            stats["sequences_stored"] += n_saved
            if stats["batches_executed"] % 10 == 0:
                print(
                    f"[NCBI {marker}] batches={stats['batches_executed']} stored_total={stats['sequences_stored']}",
                    file=sys.stderr,
                    flush=True,
                )
    return stats


def ingest_bold_coi_by_genera(
    conn: sqlite3.Connection,
    max_genera: int | None = None,
    timeout: int = 25,
) -> dict[str, int]:
    """Ingest public COI-5P sequences from BOLD Systems v4 Public API for target checklist genera."""
    name_map = resolve_species_name_Map(conn)
    existing_accessions = {
        r["accession"].split(".")[0]
        for r in conn.execute("SELECT accession FROM sequences WHERE marker='COI'").fetchall()
    }
    genera = sorted({r["genus_name"] for r in conn.execute("SELECT DISTINCT genus_name FROM taxa").fetchall()})
    done_batches = {
        r["batch_key"]
        for r in conn.execute("SELECT batch_key FROM ingestion_checkpoints WHERE source_db='BOLD'").fetchall()
    }

    session = requests.Session()
    session.headers.update({"User-Agent": "BengalFishRef/1.0 (FBD-NIB Bangladesh)"})
    stats = {"genera_queried": 0, "sequences_stored": 0}

    for genus in genera:
        bkey = f"BOLD::COI::{genus}"
        if bkey in done_batches:
            continue
        if max_genera is not None and stats["genera_queried"] >= max_genera:
            break
        try:
            resp = session.get(
                "https://v4.boldsystems.org/index.php/API_Public/sequence",
                params={"taxon": genus},
                timeout=timeout,
            )
            if resp.status_code != 200 or not resp.text.strip().startswith(">"):
                stats["genera_queried"] += 1
                continue

            n_saved = 0
            curr_hdr = ""
            curr_seq: list[str] = []
            entries: list[tuple[str, str]] = []
            for line in resp.text.splitlines():
                if line.startswith(">"):
                    if curr_hdr and curr_seq:
                        entries.append((curr_hdr, "".join(curr_seq)))
                    curr_hdr = line[1:].strip()
                    curr_seq = []
                else:
                    curr_seq.append(line.strip())
            if curr_hdr and curr_seq:
                entries.append((curr_hdr, "".join(curr_seq)))

            with conn:
                for hdr, seq in entries:
                    # BOLD FASTA format: >PROCESSID|Taxon Name|Marker|GenBankAcc
                    fields = [f.strip() for f in hdr.split("|")]
                    if len(fields) < 3:
                        continue
                    proc_id, rep_name, marker_code = fields[0], fields[1], fields[2]
                    gb_acc = fields[3] if len(fields) >= 4 else ""
                    if "COI" not in marker_code.upper():
                        continue
                    if gb_acc and gb_acc.split(".")[0] in existing_accessions:
                        continue
                    norm_rep = normalize_binomial(rep_name)
                    tx_meta = name_map.get(norm_rep.lower())
                    if tx_meta is None:
                        continue
                    if insert_ingested_sequence(
                        conn,
                        source_db="BOLD",
                        accession=proc_id,
                        reported_taxon=norm_rep,
                        taxon_meta=tx_meta,
                        marker="COI",
                        raw_seq=seq.replace("-", ""),
                        country_origin="BOLD_Public",
                    ):
                        n_saved += 1
                    if n_saved >= 25:
                        # Cap redundant copies per genus from BOLD to preserve lean DB footprint
                        break

                conn.execute(
                    """
                    INSERT OR REPLACE INTO ingestion_checkpoints (batch_key, source_db, marker, completed_at, records_fetched)
                    VALUES (?, 'BOLD', 'COI', ?, ?)
                    """,
                    (bkey, datetime.now(timezone.utc).isoformat(), n_saved),
                )
            stats["genera_queried"] += 1
            stats["sequences_stored"] += n_saved
            time.sleep(0.35)
        except Exception:
            stats["genera_queried"] += 1
            continue

    return stats
