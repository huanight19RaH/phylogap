#!/usr/bin/env python3
"""
test_b10k_manifest.py - QA Verification Suite for B10K Avian Manifest Data Integrity.

Validates schema, accession syntax, physical consistency, and evolutionary tiers.

Author: PhyloGap Research Team (QA & Verification Reviewer)
Date: 2026-10-23
"""

import unittest
import csv
import re
from pathlib import Path

class TestAvianB10KManifest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.manifest_path = Path("data/manifests/avian_b10k_manifest.csv")
        cls.assertTrue(cls.manifest_path.exists(), f"Manifest file missing: {cls.manifest_path}")
        
        with open(cls.manifest_path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            cls.rows = list(reader)
            cls.fieldnames = reader.fieldnames

    def test_file_not_empty(self):
        """Manifest must contain at least 7 curated entries."""
        self.assertGreaterEqual(len(self.rows), 7, "Manifest has too few records.")

    def test_required_columns_exist(self):
        """Verify presence of all mandatory metadata fields."""
        required_cols = [
            "tier", "role", "species", "common_name", "divergence_mya",
            "accession", "assembly_name", "assembly_level", "bioproject",
            "biosample", "total_sequence_length", "ungapped_length",
            "total_gap_length_bp", "number_of_gaps", "contig_n50",
            "scaffold_n50", "gc_percent", "sequencing_tech"
        ]
        for col in required_cols:
            self.assertIn(col, self.fieldnames, f"Missing required column: {col}")

    def test_ncbi_accession_format(self):
        """All accessions must follow NCBI GCF_ / GCA_ standard formatting."""
        acc_pattern = re.compile(r"^GC[FA]_\d{9}\.\d+$")
        for r in self.rows:
            acc = r["accession"]
            self.assertTrue(acc_pattern.match(acc), f"Invalid accession pattern: {acc}")

    def test_sequence_length_consistency(self):
        """Physical sequence lengths must be positive and total >= ungapped."""
        for r in self.rows:
            total_len = int(r["total_sequence_length"])
            ungapped_len = int(r["ungapped_length"])
            gap_len = int(r["total_gap_length_bp"])
            
            self.assertGreater(total_len, 500_000_000, f"Avian genome total length abnormally small: {total_len}")
            self.assertGreater(ungapped_len, 500_000_000, f"Ungapped length abnormally small: {ungapped_len}")
            self.assertGreaterEqual(total_len, ungapped_len, "Total sequence length < ungapped length")
            self.assertEqual(total_len - ungapped_len, gap_len, "Gap length != total - ungapped")

    def test_n50_and_gaps_validity(self):
        """Contig N50 must be <= Scaffold N50 and gap count >= 0."""
        for r in self.rows:
            c_n50 = int(r["contig_n50"])
            s_n50 = int(r["scaffold_n50"])
            gaps = int(r["number_of_gaps"])
            
            self.assertGreater(c_n50, 0, f"Contig N50 must be positive in {r['accession']}")
            self.assertGreater(s_n50, 0, f"Scaffold N50 must be positive in {r['accession']}")
            self.assertLessEqual(c_n50, s_n50, f"Contig N50 ({c_n50}) cannot exceed Scaffold N50 ({s_n50})")
            self.assertGreaterEqual(gaps, 0, f"Gap count cannot be negative in {r['accession']}")

    def test_gc_percentage_range(self):
        """Vertebrate avian GC% should fall realistically between 35% and 55%."""
        for r in self.rows:
            gc = float(r["gc_percent"])
            self.assertGreaterEqual(gc, 35.0, f"GC% abnormally low ({gc}) for {r['species']}")
            self.assertLessEqual(gc, 55.0, f"GC% abnormally high ({gc}) for {r['species']}")

    def test_evolutionary_tiers_coverage(self):
        """Ensure all 4 evolutionary tiers (RQ3) are represented."""
        tiers = {r["tier"] for r in self.rows}
        self.assertTrue(any("Tier 1" in t for t in tiers), "Tier 1 (Sister species) missing")
        self.assertTrue(any("Tier 2" in t for t in tiers), "Tier 2 (Same family) missing")
        self.assertTrue(any("Tier 3" in t for t in tiers), "Tier 3 (Same order) missing")
        self.assertTrue(any("Tier 4" in t for t in tiers), "Tier 4 (Distant order) missing")

    def test_divergence_time_ordering(self):
        """Validate evolutionary divergence progression across reference species."""
        tier1_refs = [float(r["divergence_mya"]) for r in self.rows if "Tier 1" in r["tier"] and r["role"] != "Target Species"]
        tier2_refs = [float(r["divergence_mya"]) for r in self.rows if "Tier 2:" in r["tier"] and r["role"] != "Target Species"]
        tier3_refs = [float(r["divergence_mya"]) for r in self.rows if "Tier 3" in r["tier"] and r["role"] != "Target Species"]
        tier4_refs = [float(r["divergence_mya"]) for r in self.rows if "Tier 4" in r["tier"] and r["role"] != "Target Species"]
        
        self.assertLess(tier1_refs[0], tier2_refs[0], "Tier 1 divergence must be < Tier 2")
        self.assertLess(tier2_refs[0], tier3_refs[0], "Tier 2 divergence must be < Tier 3")
        self.assertLess(tier3_refs[0], tier4_refs[0], "Tier 3 divergence must be < Tier 4")

if __name__ == "__main__":
    unittest.main()
