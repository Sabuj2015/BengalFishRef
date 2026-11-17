# BengalFishRef (`v1.0.0`)

**A Phylogenetically Audited Mitochondrial Reference Atlas, Metabarcoding Primer Benchmark, and Serverless Identification Platform for Fishes and Commercial Crustaceans of Bangladesh and the Bay of Bengal**

Developed by **Sabuj Biswas**, Fisheries Biotechnology Division, National Institute of Biotechnology (FBD-NIB), Bangladesh.

---

## Key Highlights (`v1.0.0`)
- **859 Harmonized Valid Species:** `734` Actinopterygii, `79` Elasmobranchii, and `46` commercial Decapoda across the Ganges–Brahmaputra–Meghna (GBM) delta, Sundarbans mangrove estuary, and northern Bay of Bengal (`57` historical South Asian synonyms resolved).
- **27,467 Audited Mitochondrial Accessions (`COI`, `12S`, `16S`, `Cytb`, `Mitogenomes`):** Two-stage hygiene, genetic-code-aware NUMT screening (`transl_table=2` vs. `5`), and `MAFFT` + K2P phylogenetic concordance purged **1,032 defective/mislabeled public accessions** (`222` internal stop-codon NUMTs and `638` phylogenetically mislabeled GenBank vouchers), yielding **26,435 Gold/Silver production sequences**.
- **14-Primer *In Silico* eDNA Benchmark (`25,300` Tests):** Quantifies 3'-weighted mismatch penalties, species-level amplicon collisions, and Bengal-adapted degenerate primers (`Teleo-Bengal`, `MiFish-U-Bengal`).
- **Serverless Browser Portal (`./web`):** Zero-backend species/BDTPI Dark-Taxa explorer and client-side FASTA/Sanger `.ab1` identifier + mitochondrial NUMT quality-control checker.

## Quick Start (`uv`)

```bash
# Install isolated reproducible environment
uv sync

# Run deterministic unit tests and 6-species end-to-end smoke test
uv run pytest -v
uv run bengalfishref smoke-test --threads 4

# Validate any new wet-lab FASTA or Sanger .ab1 chromatogram file
uv run bengalfishref validate sample_barcodes.fasta --output sample_qc.json

# Export QIIME2, DADA2, SINTAX, NCBI BLAST+, and static web portal
uv run bengalfishref export-all
```

## Deploying the Serverless Web Portal Free on GitHub Pages
1. Create a new public GitHub repository named `BengalFishRef` and push this folder (`main` branch).
2. In GitHub **Settings → Pages → Build and deployment**, set **Source** to **GitHub Actions**.
3. `.github/workflows/deploy-pages.yml` automatically runs `uv run pytest -v` and publishes `./web` to `https://sabuj2015.github.io/BengalFishRef/`.

## Citation & License
Code released under the **MIT License**; curated reference metadata and taxonomy distributions released under **CC-BY 4.0**.
