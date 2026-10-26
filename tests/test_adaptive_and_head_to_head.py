"""
Unit and Integration Tests for PhyloGap Phase 5:
1. Adaptive Heterogeneous Evolutionary Sequence Decay Model.
2. Head-to-Head Comparative Paradigms (DENTIST, GapPredict, Raw Homology, PhyloGap).
"""

import pytest
import random
from scripts.adaptive_evolutionary_model import (
    AdaptiveEvolutionaryDecayModel,
    BASE_LAMBDA_0,
    BASE_T_CRIT
)
from scripts.benchmark_head_to_head_sota import (
    BaselineToolSimulator,
    generate_synthetic_sparse_reads
)


class TestAdaptiveEvolutionaryModel:
    """Validates fine-grained locus-specific evolutionary scaling."""

    def setup_method(self):
        self.model = AdaptiveEvolutionaryDecayModel(base_lambda=BASE_LAMBDA_0)

    def test_feature_rate_hierarchy(self):
        """Purifying selection must maintain: exon < utr < intron < intergenic < repeat."""
        lambda_exon = self.model.calculate_effective_lambda("coding_exon")
        lambda_utr = self.model.calculate_effective_lambda("utr")
        lambda_intron = self.model.calculate_effective_lambda("intron")
        lambda_intergenic = self.model.calculate_effective_lambda("intergenic")
        lambda_repeat = self.model.calculate_effective_lambda("repeat")

        assert lambda_exon < lambda_utr
        assert lambda_utr < lambda_intron
        assert lambda_intron < lambda_intergenic
        assert lambda_intergenic < lambda_repeat

    def test_gc_content_penalty(self):
        """Extreme GC bias must elevate effective mutation rate."""
        f_neutral = self.model.compute_gc_factor(0.42)
        f_high_gc = self.model.compute_gc_factor(0.70)
        f_low_gc = self.model.compute_gc_factor(0.15)

        assert f_neutral == pytest.approx(1.0, abs=1e-4)
        assert f_high_gc > 1.10
        assert f_low_gc > 1.10

    def test_adaptive_t_crit_ordering(self):
        """Exons should remain reconstructible longer (higher t_crit) than neutral introns."""
        t_crit_exon = self.model.calculate_effective_t_crit("coding_exon")
        t_crit_intron = self.model.calculate_effective_t_crit("intron")
        t_crit_repeat = self.model.calculate_effective_t_crit("repeat")

        assert t_crit_exon > t_crit_intron
        assert t_crit_intron > t_crit_repeat
        assert t_crit_exon > 20.0  # Exons tolerate >20 Mya divergence

    def test_adaptive_homology_weight_decay(self):
        """Weight must be positive below t_crit and strictly zero at or above t_crit."""
        w_near = self.model.compute_adaptive_homology_weight(5.0, "intron")
        w_exceeded = self.model.compute_adaptive_homology_weight(50.0, "intron")

        assert w_near > 0.0
        assert w_exceeded == 0.0

    def test_classify_locus_risk_dict(self):
        """Risk classification returns valid schema and correct risk tier."""
        risk_info = self.model.classify_locus_risk(5.2, "coding_exon", 0.42)
        assert risk_info["risk_tier"] == "SAFE_CONGENERIC"
        assert "lambda_eff" in risk_info
        assert "effective_t_crit" in risk_info
        assert not risk_info["is_exceeding_t_crit"]


class TestHeadToHeadSimulators:
    """Validates baseline simulation invariants."""

    def test_dentist_zero_reads_refusal(self):
        """DENTIST must refuse to close gaps when no reads are provided."""
        rng = random.Random(42)
        res = BaselineToolSimulator.run_dentist_paradigm("ACGT" * 10, "zero_reads", rng)
        assert not res["is_closed"]
        assert "N" in res["sequence"]

    def test_gappredict_execution(self):
        """GapPredict produces valid sequence matching expected length."""
        rng = random.Random(42)
        res = BaselineToolSimulator.run_gappredict_paradigm("ACGT" * 15, "A" * 100, "C" * 100, rng)
        assert res["is_closed"]
        assert len(res["sequence"]) == 60

    def test_phylogap_rejection_beyond_tcrit(self):
        """PhyloGap must reject distant orders (>46 Mya) without reads via selective gating."""
        flank_l = "ACGTACGTAC" * 10
        flank_r = "TGCATGCATG" * 10
        gt_gap = "AAAAAAAAAA" * 10
        ref_scaffold = flank_l + ("GGGGGGGGGG" * 10) + flank_r

        res = BaselineToolSimulator.run_phylogap_paradigm(
            flank_left=flank_l,
            flank_right=flank_r,
            homology_reference=ref_scaffold,
            divergence_mya=46.2,
            expected_len=100,
            sparse_reads=[],
            ground_truth=gt_gap
        )
        assert not res["is_closed"]
        assert res["status"] == "REJECTED_SELECTIVE_GATING"
        assert res["sequence"] == "N" * 100
