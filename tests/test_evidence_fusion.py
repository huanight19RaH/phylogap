"""
Unit test suite for PhyloGap Evidence-Weighted Multimodal Fusion and Selective Gating Engine (RQ4 & RQ5).
Validates dynamic homology weighting, sparse-read polishing, composite scoring, and error gating.
"""

import unittest
import os
import sys

# Ensure repository root is on sys.path
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from scripts.evidence_fusion_engine import (
    EvidenceWeightedFuser,
    SelectivePredictionGater,
    EMPIRICAL_DECAY_LAMBDA,
    EMPIRICAL_T_CRIT
)

class TestEvidenceFusionEngine(unittest.TestCase):
    def setUp(self):
        self.fuser = EvidenceWeightedFuser()
        self.gater = SelectivePredictionGater(min_confidence_threshold=0.65)

    def test_homology_weight_decay(self):
        # At t=0, weight should be 1.0
        w0 = self.fuser.compute_homology_weight(0.0)
        self.assertAlmostEqual(w0, 1.0, places=2)

        # At t=5.2 Mya (Tier 1), weight should be moderate (~0.5 - 0.7)
        w5 = self.fuser.compute_homology_weight(5.2)
        self.assertTrue(0.40 <= w5 <= 0.70)

        # At t=11.5 Mya (Tier 2), weight should be low (~0.1 - 0.25)
        w11 = self.fuser.compute_homology_weight(11.5)
        self.assertTrue(0.05 <= w11 <= 0.25)

        # At t >= t_crit (13.2 Mya), weight should strictly drop to 0.0
        w_crit = self.fuser.compute_homology_weight(13.2)
        w_distant = self.fuser.compute_homology_weight(46.2)
        self.assertEqual(w_crit, 0.0)
        self.assertEqual(w_distant, 0.0)

    def test_sparse_read_polishing_corrects_substitution(self):
        # Homology has an evolutionary substitution: 'A' -> 'G' at index 4
        # Target ground truth is 'ACGTATGC'
        # Homology candidate is 'ACGTGTGC'
        homology_candidate = "ACGTGTGC"
        # Two sparse reads covering the locus with the true target base 'A'
        sparse_reads = ["ACGTATGC", "CGTATGC"]
        
        polished, meta = self.fuser.polish_homology_with_sparse_reads(
            homology_seq=homology_candidate,
            reads=sparse_reads,
            divergence_mya=11.5 # Homology weight is small (~0.12), read vote (2.0) easily overrides
        )

        self.assertEqual(polished, "ACGTATGC")
        self.assertEqual(meta["substitutions_corrected"], 1)
        self.assertGreater(meta["read_support_rate"], 0.70)

    def test_fusion_composite_confidence(self):
        homology_candidate = "ACGTACGTACGT"
        sparse_reads = ["ACGTACGTACGT", "ACGTACGT"]
        
        res = self.fuser.fuse(
            homology_seq=homology_candidate,
            reads=sparse_reads,
            divergence_mya=5.2,
            flank_identity=0.96,
            expected_gap_len=12
        )

        self.assertEqual(res["candidate_sequence"], homology_candidate)
        self.assertGreater(res["confidence_score"], 0.75)
        self.assertEqual(res["length_factor"], 1.0)

    def test_selective_gating_accepts_good_prediction(self):
        fusion_res = {
            "candidate_sequence": "ACGTACGTACGT",
            "confidence_score": 0.88
        }
        gate_res = self.gater.evaluate_gating(
            fusion_result=fusion_res,
            divergence_mya=5.2,
            expected_gap_len=12,
            reads_present=True
        )
        self.assertTrue(gate_res["is_accepted"])
        self.assertEqual(gate_res["final_sequence"], "ACGTACGTACGT")
        self.assertIsNone(gate_res["rejection_reason"])

    def test_selective_gating_rejects_distant_species_without_reads(self):
        # Divergence is 46.2 Mya (> t_crit=13.2 Mya) and reads_present is False
        fusion_res = {
            "candidate_sequence": "ACGTACGTACGT",
            "confidence_score": 0.50
        }
        gate_res = self.gater.evaluate_gating(
            fusion_result=fusion_res,
            divergence_mya=46.2,
            expected_gap_len=12,
            reads_present=False
        )
        self.assertFalse(gate_res["is_accepted"])
        self.assertEqual(gate_res["final_sequence"], "N" * 12)
        self.assertIn("exceeds critical threshold", gate_res["rejection_reason"])

    def test_selective_gating_rejects_excessive_length_discordance(self):
        fusion_res = {
            "candidate_sequence": "ACGTACGT" * 4, # 32 bp vs expected 10 bp (220% error)
            "confidence_score": 0.80
        }
        gate_res = self.gater.evaluate_gating(
            fusion_result=fusion_res,
            divergence_mya=5.2,
            expected_gap_len=10,
            reads_present=True
        )
        self.assertFalse(gate_res["is_accepted"])
        self.assertEqual(gate_res["final_sequence"], "N" * 10)
        self.assertIn("Length error exceeds", gate_res["rejection_reason"])

if __name__ == "__main__":
    unittest.main()
