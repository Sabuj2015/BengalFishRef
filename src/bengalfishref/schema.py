"""SQLite data layer and Bengal Dark-Taxa Priority Index (BDTPI) schema for BengalFishRef."""

from __future__ import annotations

import json
import sqlite3
from pathlib import Path
from typing import Any, Iterable

IUCN_WEIGHTS: dict[str, float] = {
    "CR": 1.00,
    "EN": 0.85,
    "VU": 0.70,
    "NT": 0.50,
    "DD": 0.45,
    "NE": 0.30,
    "LC": 0.20,
}

COMMERCIAL_WEIGHTS: dict[str, float] = {
    "Core_Commercial": 1.00,
    "Artisanal_Fishery": 0.65,
    "Conservation_Endemic": 0.80,
    "Bycatch_NonCommercial": 0.25,
}

CORE_MARKERS: tuple[str, ...] = ("COI", "12S", "16S", "Cytb", "Mitogenome")


SCHEMA_SQL = """
PRAGMA journal_mode = WAL;
PRAGMA foreign_keys = ON;

CREATE TABLE IF NOT EXISTS taxa (
    taxon_id INTEGER PRIMARY KEY,
    valid_name TEXT NOT NULL UNIQUE,
    authority TEXT DEFAULT '',
    common_name_en TEXT DEFAULT '',
    common_name_bn TEXT DEFAULT '',
    tax_group TEXT NOT NULL,
    order_name TEXT NOT NULL,
    family_name TEXT NOT NULL,
    genus_name TEXT NOT NULL,
    synonyms_json TEXT NOT NULL DEFAULT '[]',
    habitat_zones TEXT NOT NULL,
    aquaculture_tier TEXT NOT NULL,
    iucn_bd_status TEXT NOT NULL DEFAULT 'NE',
    genetic_code INTEGER NOT NULL CHECK (genetic_code IN (2, 5))
);

CREATE TABLE IF NOT EXISTS sequences (
    seq_uid TEXT PRIMARY KEY,
    source_db TEXT NOT NULL,
    accession TEXT NOT NULL,
    reported_taxon TEXT NOT NULL,
    resolved_taxon_id INTEGER REFERENCES taxa(taxon_id),
    marker TEXT NOT NULL,
    seq_length INTEGER NOT NULL,
    ambig_pct REAL NOT NULL,
    gc_pct REAL NOT NULL DEFAULT 0.0,
    country_origin TEXT DEFAULT 'Unknown',
    is_bangladesh_voucher INTEGER NOT NULL DEFAULT 0,
    has_stop_codon_numt INTEGER NOT NULL DEFAULT 0,
    homopolymer_max INTEGER NOT NULL DEFAULT 0,
    k2p_nearest_conspecific REAL,
    k2p_nearest_heterospecific REAL,
    nearest_heterospecific_taxon TEXT,
    phylo_concordance_flag TEXT DEFAULT 'Unchecked',
    quality_tier TEXT NOT NULL CHECK (quality_tier IN ('Gold', 'Silver', 'Quarantined')),
    quarantine_reason TEXT DEFAULT '',
    sequence TEXT NOT NULL
);

CREATE INDEX IF NOT EXISTS idx_seq_taxon_marker ON sequences(resolved_taxon_id, marker, quality_tier);
CREATE INDEX IF NOT EXISTS idx_seq_marker_tier ON sequences(marker, quality_tier);

CREATE TABLE IF NOT EXISTS ingestion_checkpoints (
    batch_key TEXT PRIMARY KEY,
    source_db TEXT NOT NULL,
    marker TEXT NOT NULL,
    completed_at TEXT NOT NULL,
    records_fetched INTEGER NOT NULL DEFAULT 0
);

CREATE TABLE IF NOT EXISTS primer_evaluations (
    eval_id INTEGER PRIMARY KEY AUTOINCREMENT,
    primer_name TEXT NOT NULL,
    marker TEXT NOT NULL,
    taxon_id INTEGER NOT NULL REFERENCES taxa(taxon_id),
    seq_uid TEXT NOT NULL REFERENCES sequences(seq_uid),
    fwd_mismatches INTEGER NOT NULL,
    rev_mismatches INTEGER NOT NULL,
    fwd_3prime_mm INTEGER NOT NULL,
    rev_3prime_mm INTEGER NOT NULL,
    weighted_penalty REAL NOT NULL,
    delta_tm_est REAL NOT NULL,
    amplifies_in_silico INTEGER NOT NULL,
    amplicon_length INTEGER,
    amplicon_seq TEXT,
    UNIQUE(primer_name, seq_uid)
);

CREATE TABLE IF NOT EXISTS benchmark_results (
    run_id TEXT NOT NULL,
    ecosystem TEXT NOT NULL,
    marker TEXT NOT NULL,
    error_rate REAL NOT NULL,
    db_condition TEXT NOT NULL,
    classifier TEXT NOT NULL,
    n_queries INTEGER NOT NULL,
    precision_macro REAL NOT NULL,
    recall_macro REAL NOT NULL,
    f1_macro REAL NOT NULL,
    f1_ci_low REAL NOT NULL,
    f1_ci_high REAL NOT NULL,
    false_positive_exotic_rate REAL NOT NULL,
    latency_ms_per_query REAL NOT NULL,
    PRIMARY KEY (run_id, ecosystem, marker, error_rate, db_condition, classifier)
);
"""


def connect_db(db_path: str | Path) -> sqlite3.Connection:
    """Open a WAL-enabled SQLite connection and ensure schema exists."""
    path = Path(db_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(str(path))
    conn.row_factory = sqlite3.Row
    conn.executescript(SCHEMA_SQL)
    return conn


def upsert_taxa(conn: sqlite3.Connection, rows: Iterable[dict[str, Any]]) -> int:
    """Insert or update harmonized taxonomic checklist records."""
    count = 0
    with conn:
        for row in rows:
            synonyms = row.get("synonyms", [])
            if isinstance(synonyms, list):
                synonyms_json = json.dumps(sorted(set(synonyms)))
            else:
                synonyms_json = str(synonyms)
            conn.execute(
                """
                INSERT INTO taxa (
                    taxon_id, valid_name, authority, common_name_en, common_name_bn,
                    tax_group, order_name, family_name, genus_name, synonyms_json,
                    habitat_zones, aquaculture_tier, iucn_bd_status, genetic_code
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                ON CONFLICT(valid_name) DO UPDATE SET
                    taxon_id = excluded.taxon_id,
                    authority = CASE WHEN excluded.authority != '' THEN excluded.authority ELSE taxa.authority END,
                    common_name_en = CASE WHEN excluded.common_name_en != '' THEN excluded.common_name_en ELSE taxa.common_name_en END,
                    common_name_bn = CASE WHEN excluded.common_name_bn != '' THEN excluded.common_name_bn ELSE taxa.common_name_bn END,
                    tax_group = excluded.tax_group,
                    order_name = excluded.order_name,
                    family_name = excluded.family_name,
                    genus_name = excluded.genus_name,
                    synonyms_json = excluded.synonyms_json,
                    habitat_zones = excluded.habitat_zones,
                    aquaculture_tier = excluded.aquaculture_tier,
                    iucn_bd_status = excluded.iucn_bd_status,
                    genetic_code = excluded.genetic_code
                """,
                (
                    int(row["taxon_id"]),
                    row["valid_name"].strip(),
                    row.get("authority", ""),
                    row.get("common_name_en", ""),
                    row.get("common_name_bn", ""),
                    row["tax_group"],
                    row["order_name"],
                    row["family_name"],
                    row["genus_name"],
                    synonyms_json,
                    row["habitat_zones"],
                    row["aquaculture_tier"],
                    row.get("iucn_bd_status", "NE"),
                    int(row["genetic_code"]),
                ),
            )
            count += 1
    return count


def compute_bdtpi(
    iucn_status: str,
    aquaculture_tier: str,
    gold_marker_counts: dict[str, int],
) -> float:
    """Compute the Bengal Dark-Taxa Priority Index (BDTPI in [0, 1]).

    BDTPI = 0.40 * S_conservation + 0.35 * S_commercial + 0.25 * S_marker_deficit
    Higher scores prioritize threatened and economically vital Bangladeshi taxa
    lacking Gold-tier reference coverage across core barcoding/eDNA loci.
    """
    s_cons = IUCN_WEIGHTS.get(iucn_status.upper(), 0.30)
    s_comm = COMMERCIAL_WEIGHTS.get(aquaculture_tier, 0.25)
    # Weight core eDNA/barcoding markers: COI (0.30), 12S (0.30), 16S (0.20), Mitogenome (0.20)
    marker_weights = {"COI": 0.30, "12S": 0.30, "16S": 0.20, "Mitogenome": 0.20}
    deficit = 0.0
    for m, w in marker_weights.items():
        n_gold = gold_marker_counts.get(m, 0)
        if n_gold == 0:
            deficit += w
        elif n_gold < 3 and m != "Mitogenome":
            deficit += 0.4 * w
    score = 0.40 * s_cons + 0.35 * s_comm + 0.25 * deficit
    return round(min(max(score, 0.0), 1.0), 4)
