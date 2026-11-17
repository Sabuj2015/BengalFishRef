"""Zero-leak Zenodo Deposit & DOI Publisher for BengalFishRef v1.0."""

from __future__ import annotations

import json
from pathlib import Path
import sys
from typing import Any

import dotenv
import requests

ENV_PATH = Path("/home/sabuj/.env")
WORKSPACE = Path("/home/sabuj/Research/01_Active_Projects/Antigravity")
ZIP_PATH = WORKSPACE / "02_Data" / "processed" / "releases" / "BengalFishRef_v1.0_Zenodo_Deposit.zip"
ZENODO_API = "https://zenodo.org/api/deposit/depositions"


def _load_zenodo_token() -> str:
    vals = dotenv.dotenv_values(str(ENV_PATH))
    tok = (vals.get("ZENODO_ACCESS_TOKEN") or "").strip()
    if not tok:
        raise RuntimeError("ZENODO_ACCESS_TOKEN is empty or missing in /home/sabuj/.env")
    return tok


def upload_and_publish_zenodo(zip_path: Path = ZIP_PATH, publish: bool = True) -> dict[str, Any]:
    token = _load_zenodo_token()
    if not zip_path.exists():
        raise FileNotFoundError(f"Archive not found: {zip_path}")

    session = requests.Session()
    session.headers.update(
        {
            "Authorization": f"Bearer {token}",
            "User-Agent": "BengalFishRef/1.0 (FBD-NIB Bangladesh)",
        }
    )

    # 1. Reuse existing unsubmitted BengalFishRef draft if one exists; else create a fresh deposition
    r_list = session.get(ZENODO_API, params={"status": "draft"}, timeout=30)
    if r_list.status_code != 200:
        raise RuntimeError(f"Zenodo auth/list failed ({r_list.status_code}): {r_list.text[:300]}")

    dep: dict[str, Any] | None = None
    for item in r_list.json():
        title_existing = (item.get("metadata") or {}).get("title", "")
        if "BengalFishRef" in title_existing or not title_existing:
            dep = item
            break

    if dep is None:
        r_create = session.post(
            ZENODO_API,
            headers={"Content-Type": "application/json"},
            data=json.dumps({}),
            timeout=30,
        )
        if r_create.status_code not in (200, 201):
            raise RuntimeError(f"Zenodo create deposition failed ({r_create.status_code}): {r_create.text[:300]}")
        dep = r_create.json()

    dep_id = dep["id"]
    bucket_url = dep["links"]["bucket"]

    # 2. Upload zip file to bucket if not already uploaded with identical size
    existing_files = {f.get("filename"): f.get("filesize") for f in dep.get("files", [])}
    target_size = zip_path.stat().st_size
    if existing_files.get(zip_path.name) != target_size:
        print(f"Uploading {zip_path.name} ({target_size / (1024*1024):.2f} MB) to Zenodo bucket...", file=sys.stderr, flush=True)
        with zip_path.open("rb") as fh:
            r_up = session.put(
                f"{bucket_url}/{zip_path.name}",
                data=fh,
                headers={"Content-Type": "application/octet-stream"},
                timeout=600,
            )
        if r_up.status_code not in (200, 201):
            raise RuntimeError(f"Zenodo bucket upload failed ({r_up.status_code}): {r_up.text[:300]}")

    # 3. Attach full scholarly metadata
    metadata_payload = {
        "metadata": {
            "title": "BengalFishRef v1.0: A Phylogenetically Audited Mitochondrial Reference Atlas and eDNA Benchmark for Fishes and Commercial Crustaceans of Bangladesh and the Bay of Bengal",
            "upload_type": "dataset",
            "description": (
                "<p><strong>BengalFishRef v1.0</strong> is an open-access, phylogenetically audited mitochondrial "
                "reference atlas (<code>COI</code>, <code>12S rRNA</code>, <code>16S rRNA</code>, <code>Cytb</code>, "
                "and complete mitogenomes) covering <strong>859 documented fish and commercial crustacean species</strong> "
                "of Bangladesh, the Sundarbans mangrove estuary, and the northern Bay of Bengal, developed at the "
                "<strong>Fisheries Biotechnology Division, National Institute of Biotechnology (NIB), Bangladesh</strong>.</p>"
                "<p>This release archive contains <strong>26,435 Gold/Silver curated sequences</strong> formatted for "
                "<strong>QIIME2</strong>, <strong>DADA2</strong>, <strong>SINTAX</strong>, and <strong>NCBI BLAST+</strong>, "
                "the harmonized 859-species master checklist, SHA-256 integrity manifest, and Supplementary Tables S1–S4 "
                "(BDTPI Dark-Taxa rankings, two-stage NUMT and phylogenetic mislabeling audits, 14-primer thermodynamic "
                "benchmarks, and 4-ecosystem simulated eDNA classifier evaluations with 95% bootstrap confidence intervals).</p>"
                "<p><strong>Live Serverless Web Portal:</strong> <a href='https://sabuj2015.github.io/BengalFishRef/'>https://sabuj2015.github.io/BengalFishRef/</a><br/>"
                "<strong>Source Code Repository:</strong> <a href='https://github.com/Sabuj2015/BengalFishRef'>https://github.com/Sabuj2015/BengalFishRef</a></p>"
            ),
            "creators": [
                {
                    "name": "Biswas, Sabuj",
                    "affiliation": "Fisheries Biotechnology Division, National Institute of Biotechnology (NIB), Bangladesh",
                    "orcid": "0000-0003-3320-0714",
                }
            ],
            "access_right": "open",
            "license": "cc-by-4.0",
            "version": "1.0.0",
            "keywords": [
                "environmental DNA (eDNA)",
                "DNA barcoding",
                "metabarcoding",
                "Bangladesh",
                "Bay of Bengal",
                "Ganges-Brahmaputra-Meghna delta",
                "Ichthyofauna",
                "Decapoda",
                "reference database curation",
            ],
            "related_identifiers": [
                {
                    "identifier": "https://github.com/Sabuj2015/BengalFishRef",
                    "relation": "isSupplementedBy",
                    "resource_type": "software",
                },
                {
                    "identifier": "https://sabuj2015.github.io/BengalFishRef/",
                    "relation": "isDocumentedBy",
                    "resource_type": "other",
                },
            ],
        }
    }
    r_meta = session.put(
        f"{ZENODO_API}/{dep_id}",
        headers={"Content-Type": "application/json"},
        data=json.dumps(metadata_payload),
        timeout=30,
    )
    if r_meta.status_code != 200:
        raise RuntimeError(f"Zenodo metadata update failed ({r_meta.status_code}): {r_meta.text[:400]}")

    updated = r_meta.json()
    prereserve_doi = ((updated.get("metadata") or {}).get("prereserve_doi") or {}).get("doi", "")

    if not publish:
        return {
            "status": "draft_staged",
            "deposition_id": dep_id,
            "prereserved_doi": prereserve_doi,
            "draft_url": f"https://zenodo.org/uploads/{dep_id}",
        }

    # 4. Publish to mint permanent DOI!
    r_pub = session.post(f"{ZENODO_API}/{dep_id}/actions/publish", timeout=60)
    if r_pub.status_code not in (200, 202):
        raise RuntimeError(f"Zenodo publish failed ({r_pub.status_code}): {r_pub.text[:400]}")

    pub_data = r_pub.json()
    doi = pub_data.get("doi") or (pub_data.get("metadata") or {}).get("doi") or prereserve_doi
    rec_id = pub_data.get("record_id") or pub_data.get("id") or dep_id
    return {
        "status": "published",
        "record_id": rec_id,
        "doi": doi,
        "doi_url": f"https://doi.org/{doi}" if doi else "",
        "record_html": pub_data.get("links", {}).get("record_html") or f"https://zenodo.org/records/{rec_id}",
    }


if __name__ == "__main__":
    out = upload_and_publish_zenodo(publish=True)
    print(json.dumps(out, indent=2))
