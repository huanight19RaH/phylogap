"""
Unit test suite for PhyloGap Cross-Species Homology Gap Reconstruction Engine.
Validates Synteny Anchor Mapping, Strand Inversion Handling,
Evidence-Aware Selective Prediction, and Avian Divergence Simulation.
"""

import unittest
import os
import sys

# Ensure repository root is on sys.path
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from scripts.cross_species_homology_engine import (
    SyntenicAnchorMapper,
    EvidenceAwareSelectivePredictor,
    AvianEvolutionarySimulator
)

class TestCrossSpeciesHomologyEngine(unittest.TestCase):
    def setUp(self):
        self.mapper = SyntenicAnchorMapper(max_flank_mismatch=0.20)
        self.predictor = EvidenceAwareSelectivePredictor(
            min_flank_identity=0.82,
            max_allowed_len_discordance=0.35,
            min_confidence_threshold=0.70
        )
        self.simulator = AvianEvolutionarySimulator(seed=123)

    def test_reverse_complement(self):
        seq = "ATGCAGTT"
        expected = "AACTGCAT"
        self.assertEqual(SyntenicAnchorMapper.reverse_complement(seq), expected)

    def test_syntenic_anchor_mapping_forward(self):
        flank_left = "ACGTGTACTGACGTATGCAT"
        gap_ground_truth = "TTTTAAAAGGGGCCCC"
        flank_right = "GCTAGCTAGCTAGCTA"
        ref_locus = "NNNN" + flank_left + gap_ground_truth + flank_right + "NNNN"

        res = self.mapper.map_flanks_to_reference(
            flank_left=flank_left,
            flank_right=flank_right,
            reference_locus=ref_locus,
            expected_gap_len=len(gap_ground_truth)
        )

        self.assertTrue(res["is_mapped"])
        self.assertEqual(res["orientation"], "forward")
        self.assertEqual(res["extracted_gap"], gap_ground_truth)
        self.assertAlmostEqual(res["left_identity"], 1.0, places=2)
        self.assertAlmostEqual(res["right_identity"], 1.0, places=2)
        self.assertAlmostEqual(res["length_discordance"], 0.0, places=2)

    def test_syntenic_anchor_mapping_reverse_complement(self):
        flank_left = "ACGTGTACTGACGTATGCAT"
        gap_ground_truth = "TTTTAAAAGGGGCCCC"
        flank_right = "GCTAGCTAGCTAGCTA"
        # Forward synthetic segment
        forward_segment = flank_left + gap_ground_truth + flank_right
        # Invert the reference locus
        ref_locus = "NNNN" + SyntenicAnchorMapper.reverse_complement(forward_segment) + "NNNN"

        res = self.mapper.map_flanks_to_reference(
            flank_left=flank_left,
            flank_right=flank_right,
            reference_locus=ref_locus,
            expected_gap_len=len(gap_ground_truth)
        )

        self.assertTrue(res["is_mapped"])
        self.assertEqual(res["orientation"], "reverse_complement")
        self.assertEqual(res["extracted_gap"], gap_ground_truth)

    def test_syntenic_anchor_mapping_unmapped(self):
        flank_left = "ACGTGTACTGACGTATGCAT"
        flank_right = "GCTAGCTAGCTAGCTA"
        ref_locus = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"

        res = self.mapper.map_flanks_to_reference(
            flank_left=flank_left,
            flank_right=flank_right,
            reference_locus=ref_locus,
            expected_gap_len=50
        )
        self.assertFalse(res["is_mapped"])
        self.assertEqual(res["extracted_gap"], "")

    def test_selective_predictor_accepts_high_quality_candidate(self):
        mapping_result = {
            "is_mapped": True,
            "left_identity": 0.95,
            "right_identity": 0.93,
            "length_discordance": 0.05,
            "extracted_gap": "ACGTACGTACGT"
        }
        res = self.predictor.evaluate_and_filter(mapping_result, expected_gap_len=12)
        self.assertTrue(res["is_accepted"])
        self.assertEqual(res["final_prediction"], "ACGTACGTACGT")
        self.assertGreater(res["confidence_score"], 0.85)
        self.assertIsNone(res["rejection_reason"])

    def test_selective_predictor_rejects_discordant_length(self):
        mapping_result = {
            "is_mapped": True,
            "left_identity": 0.95,
            "right_identity": 0.95,
            "length_discordance": 0.60, # Far exceeds 0.35 threshold
            "extracted_gap": "ACGT" * 10
        }
        res = self.predictor.evaluate_and_filter(mapping_result, expected_gap_len=10)
        self.assertFalse(res["is_accepted"])
        self.assertEqual(res["final_prediction"], "N" * 10)
        self.assertIn("Excessive length discordance", res["rejection_reason"])

    def test_selective_predictor_rejects_divergent_flanks(self):
        mapping_result = {
            "is_mapped": True,
            "left_identity": 0.70, # Below 0.82 threshold
            "right_identity": 0.72,
            "length_discordance": 0.10,
            "extracted_gap": "ACGTACGT"
        }
        res = self.predictor.evaluate_and_filter(mapping_result, expected_gap_len=8)
        self.assertFalse(res["is_accepted"])
        self.assertEqual(res["final_prediction"], "N" * 8)
        self.assertIn("Low flank identity", res["rejection_reason"])

    def test_avian_evolutionary_divergence_monotonicity(self):
        seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT" * 10 # 480 bp
        tiers_mya = [5.2, 11.5, 46.2, 88.0]
        mutation_counts = []

        for mya in tiers_mya:
            evolved, stats = self.simulator.evolve_sequence(seq, divergence_time_mya=mya)
            total_mutations = stats["substitutions"] + stats["insertions"] + stats["deletions"]
            mutation_counts.append(total_mutations)

        # Monotonically increasing mutations with evolutionary distance
        for k in range(len(mutation_counts) - 1):
            self.assertLessEqual(mutation_counts[k], mutation_counts[k + 1])

if __name__ == "__main__":
    unittest.main()
