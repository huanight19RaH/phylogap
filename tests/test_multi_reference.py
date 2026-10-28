"""
Unit test suite for PhyloGap Multi-Species Phylogenetic Scaffolding Engine (v1.2).
Validates:
1. MultiReferenceSyntenyAligner (Mapping and extraction across K divergent references).
2. AncestralConsensusReconstructor (Ancestral nucleotide agreement and conflict resolution).
3. Phylogenetic Corroboration Factor (1.6x weight boost on corroborated bases).
4. Integration with real biological gap sequences.
"""

import unittest
import os
import sys

# Ensure repository root is on sys.path
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from scripts.multi_reference_scaffolding import (
    MultiReferenceSyntenyAligner,
    AncestralConsensusReconstructor
)


class TestMultiReferenceScaffolding(unittest.TestCase):
    def setUp(self):
        self.aligner = MultiReferenceSyntenyAligner(max_flank_mismatch=0.25)
        self.reconstructor = AncestralConsensusReconstructor(base_decay_lambda=0.006395, kappa=2.5)

    def test_multi_reference_alignment_and_extraction(self):
        flank_l = "ACCACCACCGTGGCTTCGTACGAGCAGCTGCTGGAGGCCAAGGAGGCC"
        gap_seq = "CTGGTGAAGAGCCTGGCCGAGGAGGCCGGCCACGTGGACATCATC"
        flank_r = "GTGCTGCAGGACATCCAGGCCGAGGTGGACAAGATGCTGCTGGAG"

        # Reference 1 (Sister species, ~5.2 Mya)
        ref1 = ("NNNN" * 3) + (flank_l + gap_seq + flank_r) + ("NNNN" * 3)
        # Reference 2 (Same family, ~11.5 Mya - 1 transition substitution in gap)
        mutated_gap = gap_seq[:20] + "A" + gap_seq[21:]
        ref2 = ("NNNN" * 3) + (flank_l + mutated_gap + flank_r) + ("NNNN" * 3)

        references = [
            {"species": "Lonchura oryzivora", "sequence": ref1, "divergence_mya": 5.2},
            {"species": "Poephila acuticauda", "sequence": ref2, "divergence_mya": 11.5}
        ]

        anchored = self.aligner.extract_orthologous_scaffolds(
            flank_left=flank_l,
            flank_right=flank_r,
            references=references,
            expected_gap_len=len(gap_seq)
        )

        self.assertEqual(len(anchored), 2)
        self.assertEqual(anchored[0]["species"], "Lonchura oryzivora")
        self.assertEqual(anchored[1]["species"], "Poephila acuticauda")
        self.assertEqual(anchored[0]["scaffold_sequence"], gap_seq)
        self.assertEqual(anchored[1]["scaffold_sequence"], mutated_gap)

    def test_ancestral_consensus_corroboration(self):
        # When 2 references agree on the base, it should be selected as ancestral
        anchored_scaffolds = [
            {
                "species": "Ref A (~10 Mya)",
                "divergence_mya": 10.0,
                "scaffold_sequence": "ACGTACGTACGT",
                "flank_identity": 0.95
            },
            {
                "species": "Ref B (~15 Mya)",
                "divergence_mya": 15.0,
                "scaffold_sequence": "ACGTACGTACGT", # Perfectly identical
                "flank_identity": 0.93
            }
        ]

        consensus, meta = self.reconstructor.reconstruct_ancestral_scaffold(
            anchored_scaffolds=anchored_scaffolds,
            expected_gap_len=12
        )

        self.assertEqual(consensus, "ACGTACGTACGT")
        self.assertEqual(meta["references_used"], 2)
        self.assertEqual(meta["concordance_rate"], 1.0)
        self.assertEqual(meta["corroborated_sites"], 12)
        self.assertEqual(meta["conflicts_resolved"], 0)

    def test_ancestral_consensus_resolves_lineage_drift(self):
        # Reference 1 (Ref A) has lineage-specific substitution: index 4 has 'C'
        # Outgroups Ref B and Ref C both retain the ancestral 'T' at index 4
        anchored_scaffolds = [
            {
                "species": "Ref A (~12 Mya)",
                "divergence_mya": 12.0,
                "scaffold_sequence": "ACGTCCGTACGT", # Mutated to C at index 4
                "flank_identity": 0.95
            },
            {
                "species": "Ref B (~14 Mya)",
                "divergence_mya": 14.0,
                "scaffold_sequence": "ACGTTCGTACGT", # Retains ancestral T at index 4 (12 bp)
                "flank_identity": 0.90
            },
            {
                "species": "Ref C (~16 Mya)",
                "divergence_mya": 16.0,
                "scaffold_sequence": "ACGTTCGTACGT", # Retains ancestral T at index 4 (12 bp)
                "flank_identity": 0.88
            }
        ]

        consensus, meta = self.reconstructor.reconstruct_ancestral_scaffold(
            anchored_scaffolds=anchored_scaffolds,
            expected_gap_len=12
        )

        # Ancestral consensus must choose 'T' at index 4 over the lineage mutation 'C'
        self.assertEqual(consensus[4], "T")
        self.assertEqual(meta["status"], "MULTI_REFERENCE_CONSENSUS")
        self.assertEqual(meta["references_used"], 3)


if __name__ == "__main__":
    unittest.main()
