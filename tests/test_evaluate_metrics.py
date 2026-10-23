"""
Unit tests for PhyloGap Core Evaluation Engine (scripts/evaluate_metrics.py).
Validates exact match, substitution/indel error decomposition, False Closure detection, and AU-RC.
"""

import unittest
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from scripts.evaluate_metrics import evaluate_single_gap, evaluate_batch, compute_risk_coverage_curve

class TestPhyloGapEvaluation(unittest.TestCase):

    def test_exact_match(self):
        seq = "ATGCGTACGTTAGCTAGCTAG"
        res = evaluate_single_gap(seq, seq)
        self.assertTrue(res["exact_match"])
        self.assertEqual(res["nucleotide_accuracy"], 1.0)
        self.assertEqual(res["error_rate"], 0.0)
        self.assertEqual(res["substitutions"], 0)
        self.assertEqual(res["insertions"], 0)
        self.assertEqual(res["deletions"], 0)
        self.assertFalse(res["is_false_closure"])

    def test_unfilled_gap(self):
        true_seq = "ATGCGTACGTTAGCTAGCTAG"
        # Tool left Ns or empty
        res_empty = evaluate_single_gap(true_seq, "")
        self.assertFalse(res_empty["is_filled"])
        self.assertFalse(res_empty["is_false_closure"])
        
        res_ns = evaluate_single_gap(true_seq, "NNNNNNNNNNNNNNNNNNNNN")
        self.assertFalse(res_ns["is_filled"])
        self.assertFalse(res_ns["is_false_closure"])

    def test_substitutions(self):
        true_seq = "AAAAAAAAAA" # 10 As
        pred_seq = "AAAAATAAAA" # 1 T substitution at index 5
        res = evaluate_single_gap(true_seq, pred_seq)
        self.assertTrue(res["is_filled"])
        self.assertFalse(res["exact_match"])
        self.assertEqual(res["substitutions"], 1)
        self.assertEqual(res["insertions"], 0)
        self.assertEqual(res["deletions"], 0)
        self.assertAlmostEqual(res["nucleotide_accuracy"], 0.9, places=2)
        self.assertAlmostEqual(res["error_rate"], 0.1, places=2)
        # 10% error is right on the boundary
        self.assertFalse(res["is_false_closure"])

    def test_false_closure_detection(self):
        true_seq = "ATGCGTACGTTAGCTAGCTAG"
        # Highly mutated prediction (>30% errors)
        bad_pred = "ATGCGGGGGTTAGCTAAAAAG"
        res = evaluate_single_gap(true_seq, bad_pred, max_error_tolerance=0.10)
        self.assertTrue(res["is_filled"])
        self.assertTrue(res["is_false_closure"])
        self.assertGreater(res["error_rate"], 0.10)

    def test_length_discrepancy_false_closure(self):
        true_seq = "ATGCGTACGTTAGCTAGCTAG" # length 21
        # Truncated prediction (length 10 -> >50% length error)
        trunc_pred = "ATGCGTACGT"
        res = evaluate_single_gap(true_seq, trunc_pred, max_len_tolerance=0.20)
        self.assertTrue(res["is_filled"])
        self.assertTrue(res["is_false_closure"])

    def test_batch_evaluation_and_risk_coverage(self):
        records = [
            {"true_seq": "ATGCATGC", "pred_seq": "ATGCATGC", "confidence": 0.99}, # exact
            {"true_seq": "ATGCATGC", "pred_seq": "ATGCATGT", "confidence": 0.85}, # 1 sub (12.5% err)
            {"true_seq": "ATGCATGC", "pred_seq": "CCCCGGGG", "confidence": 0.20}, # completely wrong
            {"true_seq": "ATGCATGC", "pred_seq": "NNNNNNNN", "confidence": 0.05}, # unfilled
        ]
        batch_res = evaluate_batch(records, max_error_tolerance=0.10)
        self.assertEqual(batch_res["total_gaps"], 4)
        self.assertEqual(batch_res["filled_count"], 3)
        self.assertEqual(batch_res["exact_match_count"], 1)
        self.assertGreaterEqual(batch_res["false_closure_count"], 1)
        self.assertIn("aurc", batch_res)
        self.assertGreaterEqual(batch_res["aurc"], 0.0)

if __name__ == "__main__":
    unittest.main()
