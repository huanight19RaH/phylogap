#!/usr/bin/env python3
"""
query_b10k_metadata.py - NCBI Datasets API Metadata Extractor for PhyloGap Avian Clades.

Queries NCBI Datasets v2alpha API to extract assembly metrics, BioProject, BioSample,
scaffold/contig N50, gap counts, and sequence lengths for curated avian species.

Author: PhyloGap Research Team (Data Curator Specialist)
Date: 2026-10-23
"""

import sys
import json
import urllib.request
import urllib.error
import csv
from pathlib import Path

CANDIDATE_ACCESSIONS = [
    {
        "tier": "Tier 1: Sister Species (~3-6 Mya)",
        "role": "Target Species",
        "species": "Lonchura striata",
        "common_name": "White-rumped munia / Bengalese finch",
        "clade_level": "Sister genus / Congeneric",
        "accession": "GCF_046129695.1",
        "divergence_mya": 0.0,
        "ref_citation": "VGP / B10K (2024)"
    },
    {
        "tier": "Tier 1: Sister Species (~3-6 Mya)",
        "role": "Homology Reference",
        "species": "Lonchura oryzivora",
        "common_name": "Java sparrow",
        "clade_level": "Sister species (Congeneric)",
        "accession": "GCA_056776775.1",
        "divergence_mya": 5.2,
        "ref_citation": "TimeTree v5 / Olsson et al. (2020)"
    },
    {
        "tier": "Tier 2: Same Family (~10-15 Mya)",
        "role": "Target Species",
        "species": "Taeniopygia guttata",
        "common_name": "Zebra finch",
        "clade_level": "Estrildidae anchor",
        "accession": "GCF_003957565.2",
        "divergence_mya": 0.0,
        "ref_citation": "VGP Phase 1 (Rhie et al. Nature 2021)"
    },
    {
        "tier": "Tier 2: Same Family (~10-15 Mya)",
        "role": "Homology Reference",
        "species": "Lonchura striata",
        "common_name": "White-rumped munia",
        "clade_level": "Same family (Estrildidae)",
        "accession": "GCF_046129695.1",
        "divergence_mya": 11.5,
        "ref_citation": "TimeTree v5 / Feng et al. Nature 2020 (B10K)"
    },
    {
        "tier": "Tier 2b: Same Family Alternate (~10-15 Mya)",
        "role": "Homology Reference Alternate",
        "species": "Poephila acuticauda",
        "common_name": "Long-tailed finch",
        "clade_level": "Same family (Estrildidae - Australian clade)",
        "accession": "GCA_057534655.1",
        "divergence_mya": 14.1,
        "ref_citation": "TimeTree v5 / Singhal et al. (2024)"
    },
    {
        "tier": "Tier 3: Same Order (~35-45 Mya)",
        "role": "Target Species",
        "species": "Gallus gallus",
        "common_name": "Chicken / Red junglefowl",
        "clade_level": "Galliformes anchor (Phasianidae)",
        "accession": "GCF_016699485.2",
        "divergence_mya": 0.0,
        "ref_citation": "GRCg7b / bGalGal1 (Warren et al. 2021)"
    },
    {
        "tier": "Tier 3: Same Order (~35-45 Mya)",
        "role": "Homology Reference",
        "species": "Numida meleagris",
        "common_name": "Helmeted guineafowl",
        "clade_level": "Same order (Galliformes: Numididae)",
        "accession": "GCF_002078875.1",
        "divergence_mya": 46.2,
        "ref_citation": "TimeTree v5 / Vignal et al. (2019)"
    },
    {
        "tier": "Tier 4: Distant Orders (>80 Mya)",
        "role": "Target Species",
        "species": "Calypte anna",
        "common_name": "Anna's hummingbird",
        "clade_level": "Neoaves (Apodiformes / Strisores)",
        "accession": "GCF_003957555.1",
        "divergence_mya": 0.0,
        "ref_citation": "VGP Phase 1 (Rhie et al. Nature 2021)"
    },
    {
        "tier": "Tier 4: Distant Orders (>80 Mya)",
        "role": "Homology Reference",
        "species": "Gallus gallus",
        "common_name": "Chicken",
        "clade_level": "Galloanserae (Galliformes)",
        "accession": "GCF_016699485.2",
        "divergence_mya": 88.0,
        "ref_citation": "TimeTree v5 / Jarvis et al. Science 2014; Stiller et al. 2024"
    }
]

def fetch_ncbi_assembly_report(accession: str) -> dict:
    """Fetch structured metadata from NCBI Datasets v2alpha API."""
    url = f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{accession}/dataset_report"
    req = urllib.request.Request(url, headers={"User-Agent": "PhyloGapResearch/1.0 (academic; b10k-curation)"})
    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            data = json.loads(resp.read().decode("utf-8"))
        reports = data.get("reports", [])
        if not reports:
            return {}
        return reports[0]
    except urllib.error.URLError as e:
        print(f"[ERROR] Failed to query NCBI Datasets API for {accession}: {e}", file=sys.stderr)
        return {}

def safe_int(val, default=0) -> int:
    try:
        if val is None:
            return default
        return int(val)
    except (ValueError, TypeError):
        return default

def safe_float(val, default=0.0) -> float:
    try:
        if val is None:
            return default
        return float(val)
    except (ValueError, TypeError):
        return default

def parse_assembly_metrics(report: dict) -> dict:
    """Extract standard metadata fields from NCBI Datasets report."""
    assm = report.get("assembly_info", {})
    stats = report.get("assembly_stats", {})
    org = report.get("organism", {})
    biosample = assm.get("biosample", {}) or report.get("biosample", {})
    
    total_len = safe_int(stats.get("total_sequence_length"))
    ungapped_len = safe_int(stats.get("total_ungapped_length"))
    n_contigs = safe_int(stats.get("number_of_contigs"))
    n_scaffolds = safe_int(stats.get("number_of_scaffolds"))
    
    # Gap count in scaffolds: contigs - scaffolds if not explicitly given
    explicit_gaps = stats.get("number_of_gaps")
    if explicit_gaps is not None:
        gap_count = safe_int(explicit_gaps)
    else:
        gap_count = max(0, n_contigs - n_scaffolds) if n_contigs >= n_scaffolds else 0
        
    total_gap_len = max(0, total_len - ungapped_len) if total_len and ungapped_len else 0
    
    # Extract BioProject
    bp = assm.get("bioproject_accession")
    if not bp and assm.get("bioproject_lineage"):
        bp = assm["bioproject_lineage"][0].get("bioproject_accession", "N/A")
    if not bp:
        bp = "N/A"
        
    acc_val = report.get("accession") or assm.get("assembly_accession") or "N/A"
    
    return {
        "accession": acc_val,
        "assembly_name": assm.get("assembly_name", "N/A"),
        "assembly_level": assm.get("assembly_level", "N/A"),
        "organism_name": org.get("organism_name", "N/A"),
        "tax_id": safe_int(org.get("tax_id")),
        "bioproject": bp,
        "biosample": biosample.get("accession", "N/A"),
        "sequencing_tech": assm.get("sequencing_tech", "N/A"),
        "submitter": assm.get("submitter", "N/A"),
        "total_sequence_length": total_len,
        "ungapped_length": ungapped_len,
        "total_gap_length_bp": total_gap_len,
        "number_of_gaps": gap_count,
        "contig_n50": safe_int(stats.get("contig_n50")),
        "scaffold_n50": safe_int(stats.get("scaffold_n50")),
        "gc_percent": safe_float(stats.get("gc_percent")),
        "submission_date": assm.get("release_date", "N/A")
    }

def main():
    print("[INFO] PhyloGap Avian Clade Curation Script Initializing...")
    output_dir = Path("data/manifests")
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_csv = output_dir / "avian_b10k_manifest.csv"
    
    records = []
    seen_accs = {}
    
    for item in CANDIDATE_ACCESSIONS:
        acc = item["accession"]
        if acc not in seen_accs:
            print(f"[FETCH] Querying NCBI for {acc} ({item['species']})...")
            rep = fetch_ncbi_assembly_report(acc)
            metrics = parse_assembly_metrics(rep) if rep else {}
            seen_accs[acc] = metrics
        else:
            metrics = seen_accs[acc]
            
        combined = {**item, **metrics}
        records.append(combined)
        
    print(f"[INFO] Successfully queried {len(seen_accs)} unique assemblies.")
    
    # Write to CSV
    fieldnames = [
        "tier", "role", "species", "common_name", "clade_level", "divergence_mya",
        "accession", "assembly_name", "assembly_level", "bioproject", "biosample",
        "tax_id", "total_sequence_length", "ungapped_length", "total_gap_length_bp",
        "number_of_gaps", "contig_n50", "scaffold_n50", "gc_percent", "sequencing_tech",
        "submitter", "submission_date", "ref_citation"
    ]
    
    with open(manifest_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for r in records:
            writer.writerow(r)
            
    print(f"[SUCCESS] Curated manifest saved to: {manifest_csv}")
    
    # Print summary table
    print("\n" + "="*125)
    print(f"{'Species':<22} | {'Accession':<16} | {'Level':<10} | {'Contig N50':<12} | {'Scaffold N50':<12} | {'Gaps':<8} | {'Gap bp':<10} | {'Divergence':<10}")
    print("="*125)
    for r in records:
        sp = r.get("species", "N/A")
        acc = r.get("accession", "N/A")
        lvl = str(r.get("assembly_level", "N/A"))
        c_n50 = safe_int(r.get("contig_n50"))
        s_n50 = safe_int(r.get("scaffold_n50"))
        gaps = safe_int(r.get("number_of_gaps"))
        gap_bp = safe_int(r.get("total_gap_length_bp"))
        div = f"{r.get('divergence_mya', 0.0):.1f} Mya"
        print(f"{sp:<22} | {acc:<16} | {lvl:<10} | {c_n50:<12,d} | {s_n50:<12,d} | {gaps:<8,d} | {gap_bp:<10,d} | {div:<10}")
    print("="*125)

if __name__ == "__main__":
    main()
