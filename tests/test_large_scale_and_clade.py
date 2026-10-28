"""
Unit tests for high throughput alignment, stepping stone mega-gap division,
haplotype heterozygosity detection, and clade-adaptive evolutionary routing.
"""

import pytest
import numpy as np
from scripts.high_throughput_aligner import (
    banded_needleman_wunsch,
    HaplotypeHeterozygosityDetector,
    SteppingStoneDivider
)
from scripts.clade_adaptive_router import (
    CladeAdaptiveRouter,
    CladeType,
    CLADE_REGISTRY
)
from scripts.evaluate_metrics import needleman_wunsch


class TestHighThroughputAligner:
    def test_banded_needleman_wunsch_score_parity(self):
        s1 = "ATGCGATCGATCGATCGATCGATCGATC"
        s2 = "ATGCGATCCATCGATCGATCGATCGATC" # 1 SNP
        a1_std, a2_std, score_std = needleman_wunsch(s1, s2)
        a1_bnd, a2_bnd, score_bnd = banded_needleman_wunsch(s1, s2, band_width=10)

        assert score_bnd == score_std
        assert len(a1_bnd) == len(a2_bnd)

    def test_heterozygosity_and_iupac_ambiguity(self):
        consensus = "ATGCGATCGATC"
        # Reads with clear heterozygous A/G at position 4 (G vs A)
        reads = [
            "ATGCGATCGATC",
            "ATGCGATCGATC",
            "ATGCCATCGATC", # secondary allele C at pos 4
            "ATGCCATCGATC",
        ]
        res = HaplotypeHeterozygosityDetector.audit_heterozygosity(
            reads, consensus, min_alt_depth=2, min_alt_freq=0.25
        )
        assert res["is_heterozygous"] is True
        assert res["heterozygous_sites_count"] >= 1
        assert "ambiguity_consensus" in res
        # G + C -> S in IUPAC
        assert 'S' in res["ambiguity_consensus"]

    def test_stepping_stone_divider(self):
        # 600 bp synthetic mega gap
        divider = SteppingStoneDivider(kmer_size=15, max_subgap_len=200)
        anchor_kmer = "CGACTAGCTAGCTAG" # 15 bp
        homology = ("A" * 150) + anchor_kmer + ("T" * 200) + anchor_kmer + ("G" * 150)
        reads = [
            ("A" * 20) + anchor_kmer + ("T" * 20),
            ("T" * 20) + anchor_kmer + ("G" * 20)
        ]
        part = divider.partition_mega_gap(gap_len=len(homology), homology_seq=homology, reads=reads)
        assert "is_partitioned" in part
        assert "subgaps" in part
        assert len(part["subgaps"]) >= 1


class TestCladeAdaptiveRouter:
    def test_clade_detection(self):
        router = CladeAdaptiveRouter()
        p_bird = router.detect_clade_from_name_or_taxid("Calypte anna (Anna's hummingbird)")
        p_mammal = router.detect_clade_from_name_or_taxid("Homo sapiens (Human)")
        p_plant = router.detect_clade_from_name_or_taxid("Arabidopsis thaliana")
        p_insect = router.detect_clade_from_name_or_taxid("Drosophila melanogaster")

        assert p_bird.clade_type == CladeType.AVES
        assert p_mammal.clade_type == CladeType.MAMMALIA
        assert p_plant.clade_type == CladeType.VIRIDIPLANTAE
        assert p_insect.clade_type == CladeType.INSECTA

    def test_clade_transition_probabilities(self):
        router = CladeAdaptiveRouter()
        # Mammals have kappa = 3.8, Insects have kappa = 2.0
        p_trans_mam = router.calculate_clade_transition_probability(10.0, CladeType.MAMMALIA, is_transition=True)
        p_transv_mam = router.calculate_clade_transition_probability(10.0, CladeType.MAMMALIA, is_transition=False)

        assert p_trans_mam > p_transv_mam

    def test_cpg_hypermutation_audit(self):
        router = CladeAdaptiveRouter()
        # Reference has CG, target mutated to TG (C->T deamination)
        ref = "ATCGATCGAT"
        tgt = "ATTGATTGAT"
        audit = router.audit_cpg_context(ref, tgt, CladeType.MAMMALIA)

        assert audit["cpg_sites_evaluated"] == 2
        assert audit["cpg_deaminations_observed"] == 2
        assert audit["cpg_deamination_rate"] == 1.0
        assert audit["epigenetic_drift_discounted"] is True

    def test_clade_gating_and_wgd(self):
        router = CladeAdaptiveRouter()
        # Exon at 25 Mya in Bird vs Plant
        dec_bird = router.compute_clade_gating_decision(25.0, "exon", CladeType.AVES)
        dec_plant = router.compute_clade_gating_decision(25.0, "exon", CladeType.VIRIDIPLANTAE)

        assert dec_plant["wgd_paralogy_warning"] is True
        assert dec_bird["wgd_paralogy_warning"] is False
        assert dec_bird["is_safe_to_close"] is True # bird exon t_crit is 31.4
