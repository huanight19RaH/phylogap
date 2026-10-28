"""
Unit test suite for PhyloGap Biological Feature Extractor and Preprocessing Module (v1.1).
Validates:
1. FlankUniquenessFilter (K-mer entropy, homopolymer detection, anchor ambiguity).
2. SparseReadQualityFilter (Read cleaning, non-biological filter, localized k-mer binning).
3. GapBiologicalProfiler (GC-isochore classification).
4. KimuraEvolutionaryMatrix (K80 transition vs transversion dynamics).
5. CodonFrameExtractor (Triplet frame positions, reading frame preservation).
6. RichFeatureVectorBuilder (22-dimensional feature tensor).
7. End-to-end integration with EvidenceWeightedFuser and SelectivePredictionGater.
"""

import unittest
import os
import sys
import numpy as np

# Ensure repository root is on sys.path
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from scripts.biological_feature_extractor import (
    FlankUniquenessFilter,
    SparseReadQualityFilter,
    GapBiologicalProfiler,
    KimuraEvolutionaryMatrix,
    CodonFrameExtractor,
    RichFeatureVectorBuilder
)
from scripts.evidence_fusion_engine import (
    EvidenceWeightedFuser,
    SelectivePredictionGater
)
from scripts.phylogap_cli import reconstruct_gap_pipeline


class TestBiologicalFeatureExtractor(unittest.TestCase):
    def setUp(self):
        self.flank_filter = FlankUniquenessFilter(kmer_size=5, min_entropy_threshold=1.35)
        self.read_filter = SparseReadQualityFilter(min_read_len=15)
        self.k80 = KimuraEvolutionaryMatrix(base_rate=0.006395, kappa=2.5)

    def test_kmer_entropy_and_homopolymers(self):
        # Repetitive homopolymer tract has zero/near-zero entropy
        homopolymer_seq = "AAAAAAAAAAAAAAAAAAAAAAAA"
        h_entropy = self.flank_filter.calculate_kmer_entropy(homopolymer_seq)
        self.assertEqual(h_entropy, 0.0)

        runs = self.flank_filter.detect_homopolymers(homopolymer_seq, min_run_length=5)
        self.assertEqual(len(runs), 1)
        self.assertEqual(runs[0][0], "A")
        self.assertEqual(runs[0][2], 24)

        # Complex biological sequence has high entropy (>1.5)
        complex_seq = "ATGACCACCACCGTGGCTTCGTACGAGCAGCTGCTGGAGGCCAAGGAGGCC"
        c_entropy = self.flank_filter.calculate_kmer_entropy(complex_seq)
        self.assertGreater(c_entropy, 1.5)

    def test_flank_anchor_assessment(self):
        ambiguous_flank = "AAAAAAAAAAAAAAA"
        extended_context = "ATGACCACCACCGTGGCTTCGTACGAGCAGCTGCTGGAGGCCAAGGAGGCC"
        assessment = self.flank_filter.assess_flank_anchor(
            flank_seq=ambiguous_flank,
            extended_context=extended_context
        )
        self.assertTrue(assessment["expansion_applied"])
        self.assertEqual(assessment["final_flank"], extended_context)

    def test_sparse_read_quality_filter(self):
        raw_reads = [
            "ACGT", # Too short (<15 bp)
            "NNNNNNNNNNNNNNNNN", # Invalid bases
            "AAAAAAAAAAAAAAAAAAAA", # Homopolymer saturated
            "CTGGTGAAGAGCCTGGCCGAGGAGGCCGG", # Valid biological read overlapping gap
            "TTTTTGGGGGCCCCCAAAAA" # Random non-overlapping read
        ]
        flank_left = "ATGACCACCACCGTGGCTTCGTACGAGCAG"
        flank_right = "GTGCTGCAGGACATCCAGGCCGAGGTGGAC"
        homology_seq = "CTGGTGAAGAGCCTGGCCGAGGAGGCCGGCCACGTGGACATCATCACCGTG"

        cleaned = self.read_filter.filter_and_bin_reads(
            raw_reads=raw_reads,
            flank_left=flank_left,
            flank_right=flank_right,
            homology_seq=homology_seq,
            kmer_match_len=10
        )

        self.assertEqual(len(cleaned), 1)
        self.assertEqual(cleaned[0], "CTGGTGAAGAGCCTGGCCGAGGAGGCCGG")

    def test_gap_biological_profiler_isochores(self):
        # High GC sequence (Microchromosome hotspot)
        gc_rich = "GCGCGCGCGCGCGCGCGCGC"
        gc_frac = GapBiologicalProfiler.compute_gc_fraction(gc_rich)
        self.assertAlmostEqual(gc_frac, 1.0, places=2)
        self.assertEqual(GapBiologicalProfiler.categorize_isochore(gc_frac), "H3_microchromosome_hotspot")

        # AT-rich sequence (Macrochromosome L1)
        at_rich = "ATATATATATATATATATAT"
        at_frac = GapBiologicalProfiler.compute_gc_fraction(at_rich)
        self.assertEqual(GapBiologicalProfiler.categorize_isochore(at_frac), "L1_AT_rich")

    def test_kimura_k80_substitution_dynamics(self):
        # At t=0, identical base prob should be 1.0, transition/transversion 0.0
        p0 = self.k80.transition_probabilities(0.0)
        self.assertEqual(p0["p_ident"], 1.0)
        self.assertEqual(p0["p_ts"], 0.0)
        self.assertEqual(p0["p_tv"], 0.0)

        # At t=15.0 Mya, P_ident > P_transition > P_transversion
        p15 = self.k80.transition_probabilities(15.0)
        self.assertGreater(p15["p_ident"], p15["p_ts"])
        self.assertGreater(p15["p_ts"], p15["p_tv"])

        # Transition A -> G vs Transversion A -> C
        p_ts = self.k80.get_substitution_probability("A", "G", 15.0)
        p_tv = self.k80.get_substitution_probability("A", "C", 15.0)
        self.assertGreater(p_ts, p_tv)
        self.assertAlmostEqual(p_ts / max(1e-9, p_tv), self.k80.kappa, delta=0.5)

    def test_codon_frame_extractor(self):
        # Triplet position cycling: 1 -> 2 -> 3 -> 1 -> 2 -> 3
        positions = [CodonFrameExtractor.get_codon_position(i) for i in range(6)]
        self.assertEqual(positions, [1, 2, 3, 1, 2, 3])

        # Reading frame preservation: multiples of 3 are preserved
        self.assertTrue(CodonFrameExtractor.is_reading_frame_preserved(9))
        self.assertTrue(CodonFrameExtractor.is_reading_frame_preserved(0))
        self.assertFalse(CodonFrameExtractor.is_reading_frame_preserved(10))
        self.assertFalse(CodonFrameExtractor.is_reading_frame_preserved(1))

    def test_rich_feature_vector_builder(self):
        builder = RichFeatureVectorBuilder(self.k80)
        read_counts = {"A": 2.0, "C": 0.0, "G": 0.0, "T": 0.0}

        vec = builder.build_locus_feature_vector(
            homology_base="A",
            read_counts=read_counts,
            divergence_mya=11.5,
            feature_type="coding_exon",
            locus_idx=0,
            frame_offset=0,
            local_gc_fraction=0.45,
            homology_weight=0.85
        )

        self.assertEqual(vec.shape, (22,))
        # One-hot Homology base 'A' at index 0
        self.assertEqual(vec[0], 1.0)
        self.assertEqual(vec[1], 0.0)
        # Exon flag at index 12
        self.assertEqual(vec[12], 1.0)
        # Codon position 1 at index 16
        self.assertEqual(vec[16], 1.0)
        # GC fraction at index 19
        self.assertAlmostEqual(vec[19], 0.45, places=2)

    def test_end_to_end_pipeline_with_preprocessing_and_exon_gating(self):
        flank_l = "ATGACCACCACCGTGGCTTCGTACGAGCAGCTGCTGGAGGCCAAGGAGGCC"
        flank_r = "GTGCTGCAGGACATCCAGGCCGAGGTGGACAAGATGCTGCTGGAGGCCGAG"
        gap_truth = "CTGGTGAAGAGCCTGGCCGAGGAGGCCGGCCACGTGGACATCATCACCGTG"
        target_seq = flank_l + gap_truth + flank_r

        # Reference sequence with 30 bp gap
        ref_locus = ("NNNN" * 3) + target_seq + ("NNNN" * 3)

        res = reconstruct_gap_pipeline(
            flank_left=flank_l,
            flank_right=flank_r,
            homology_reference=ref_locus,
            divergence_mya=5.2,
            expected_gap_len=len(gap_truth),
            sparse_reads=[gap_truth[:30], gap_truth[20:]],
            ground_truth_gap=gap_truth,
            feature_type="coding_exon",
            enable_preprocessing=True
        )

        self.assertTrue(res["is_accepted"])
        self.assertIn("preprocessing", res)
        self.assertTrue(res["preprocessing"]["preprocessing_enabled"])
        self.assertEqual(res["preprocessing"]["feature_type"], "coding_exon")
        self.assertGreater(res["ground_truth_eval"]["nucleotide_accuracy"], 0.95)


if __name__ == "__main__":
    unittest.main()
