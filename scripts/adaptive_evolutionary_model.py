"""
PhyloGap - Adaptive Heterogeneous Evolutionary Sequence Decay Model.
Extends the baseline empirical decay model (lambda_0 = 0.006395 / Myr) to account for:
1. Genomic functional constraint (Exon / Coding, Intron, Intergenic / Repeat).
2. Local GC content composition and isochore stability.
3. Locus-specific dynamic critical threshold t_crit_eff and homology weight w_hom.

Author: PhyloGap Research Team
Date: 2026-10-26
"""

import math
from typing import Dict, Any, Tuple

# Base empirically calibrated decay rate for avian genomes
BASE_LAMBDA_0 = 0.006395
BASE_T_CRIT = 13.2

# Relative evolutionary rate scaling factors across genomic feature types
# Exons undergo purifying selection (s_feat < 1.0)
# Introns experience near-neutral drift (s_feat ~ 1.15)
# Intergenic/Repeat regions experience fast structural turnover (s_feat ~ 1.35)
FEATURE_RATE_FACTORS = {
    "coding_exon": 0.42,
    "utr": 0.75,
    "intron": 1.15,
    "intergenic": 1.35,
    "repeat": 1.55
}

# Neutral avian reference GC fraction (~42%)
AVIAN_NEUTRAL_GC = 0.42
GC_BIAS_SENSITIVITY = 0.65


class AdaptiveEvolutionaryDecayModel:
    """
    Computes fine-grained, locus-specific evolutionary decay rates
    and selective gating thresholds.
    """

    def __init__(self, base_lambda: float = BASE_LAMBDA_0):
        self.base_lambda = base_lambda

    def compute_gc_factor(self, gc_fraction: float) -> float:
        """
        Calculates mutation acceleration factor based on GC deviation from neutral avian equilibrium.
        GC-rich isochores in birds exhibit elevated recombination and mutation rates.
        """
        gc = max(0.0, min(1.0, float(gc_fraction)))
        deviation = abs(gc - AVIAN_NEUTRAL_GC)
        return 1.0 + (GC_BIAS_SENSITIVITY * deviation)

    def get_feature_factor(self, feature_type: str) -> float:
        """Returns the rate scaling factor for a specific genomic feature."""
        feat_clean = feature_type.lower().strip()
        return FEATURE_RATE_FACTORS.get(feat_clean, 1.0)

    def calculate_effective_lambda(
        self,
        feature_type: str = "intron",
        gc_fraction: float = AVIAN_NEUTRAL_GC
    ) -> float:
        """
        Calculates effective mutation rate:
        lambda_eff = lambda_0 * s_feat * f(GC)
        """
        s_feat = self.get_feature_factor(feature_type)
        f_gc = self.compute_gc_factor(gc_fraction)
        return self.base_lambda * s_feat * f_gc

    def calculate_effective_t_crit(
        self,
        feature_type: str = "intron",
        gc_fraction: float = AVIAN_NEUTRAL_GC,
        target_accuracy_threshold: float = 91.9
    ) -> float:
        """
        Calculates effective critical divergence time threshold (Mya)
        at which nucleotide identity drops below target_accuracy_threshold (default 91.9%).
        """
        lambda_eff = self.calculate_effective_lambda(feature_type, gc_fraction)
        if lambda_eff <= 0:
            return 100.0
        # Acc(t) = 100 * exp(-lambda_eff * t) => t = -ln(Acc / 100) / lambda_eff
        ratio = max(0.01, min(0.999, target_accuracy_threshold / 100.0))
        return -math.log(ratio) / lambda_eff

    def expected_identity(
        self,
        divergence_mya: float,
        feature_type: str = "intron",
        gc_fraction: float = AVIAN_NEUTRAL_GC
    ) -> float:
        """Computes expected sequence identity percentage: Acc(t) = exp(-lambda_eff * t) * 100."""
        t = max(0.0, float(divergence_mya))
        lambda_eff = self.calculate_effective_lambda(feature_type, gc_fraction)
        return math.exp(-lambda_eff * t) * 100.0

    def compute_adaptive_homology_weight(
        self,
        divergence_mya: float,
        feature_type: str = "intron",
        gc_fraction: float = AVIAN_NEUTRAL_GC
    ) -> float:
        """
        Computes dynamic homology weight w_hom(t, feat, GC).
        Decays exponentially to 0 as divergence exceeds locus-specific t_crit_eff.
        """
        t = max(0.0, float(divergence_mya))
        t_crit = self.calculate_effective_t_crit(feature_type, gc_fraction)

        if t >= t_crit:
            return 0.0

        lambda_eff = self.calculate_effective_lambda(feature_type, gc_fraction)
        # Normalized decay weight
        weight = math.exp(-lambda_eff * t)
        # Smooth penalty factor near critical threshold
        margin = max(0.0, (t_crit - t) / t_crit)
        penalty = math.tanh(3.0 * margin)

        return round(float(weight * penalty), 4)

    def classify_locus_risk(
        self,
        divergence_mya: float,
        feature_type: str = "intron",
        gc_fraction: float = AVIAN_NEUTRAL_GC
    ) -> Dict[str, Any]:
        """Provides full risk assessment for a specific locus."""
        t = max(0.0, float(divergence_mya))
        exp_id = self.expected_identity(t, feature_type, gc_fraction)
        lambda_eff = self.calculate_effective_lambda(feature_type, gc_fraction)
        t_crit = self.calculate_effective_t_crit(feature_type, gc_fraction)
        w_hom = self.compute_adaptive_homology_weight(t, feature_type, gc_fraction)

        if t <= 6.0:
            tier = "SAFE_CONGENERIC"
            label = "Safe Zone (Congeneric / Sister Species)"
            color = "#10B981"
        elif t <= t_crit:
            tier = "CAUTION_FAMILY"
            label = "Caution Zone (Family Level Clade)"
            color = "#F59E0B"
        else:
            tier = "DANGER_ORDER_CLASS"
            label = "Danger Zone (Distant Orders / Classes)"
            color = "#F43F5E"

        return {
            "divergence_mya": t,
            "feature_type": feature_type,
            "gc_fraction": gc_fraction,
            "lambda_eff": round(lambda_eff, 6),
            "effective_t_crit": round(t_crit, 2),
            "expected_identity": round(exp_id, 2),
            "adaptive_homology_weight": w_hom,
            "risk_tier": tier,
            "risk_label": label,
            "risk_color": color,
            "is_exceeding_t_crit": bool(t >= t_crit)
        }
