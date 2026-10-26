"""
PhyloGap - Evidence-Weighted Multimodal Fusion & Selective Gating Engine.
Implements the fusion of sparse read evidence, evolutionary-weighted homology scaffolds,
and local nucleotide priors to reconstruct assembly gaps with guaranteed low False-Closure Rates.

Addresses Research Questions RQ4 & RQ5:
- RQ4: Multi-Source Evidence Fusion (Reads + Homology + Priors)
- RQ5: Selective Prediction & Error Suppression
"""

import math
from typing import List, Dict, Any, Tuple, Optional
import numpy as np

# Global constants derived empirically from RQ3 avian evaluation
EMPIRICAL_DECAY_LAMBDA = 0.006395  # Decay rate per Myr (R^2 > 0.998)
EMPIRICAL_T_CRIT = 13.2             # Critical divergence time (Mya) where raw FCR exceeds 20%

class EvidenceWeightedFuser:
    """
    Fuses sparse read fragments with a homologous scaffold, dynamically weighting
    the homology contribution by evolutionary divergence distance (Mya).
    """
    def __init__(
        self,
        decay_lambda: float = EMPIRICAL_DECAY_LAMBDA,
        t_crit: float = EMPIRICAL_T_CRIT,
        min_base_qual: float = 0.60
    ):
        self.decay_lambda = decay_lambda
        self.t_crit = t_crit
        self.min_base_qual = min_base_qual

    def compute_homology_weight(self, divergence_mya: float) -> float:
        """
        Calculates dynamic homology confidence weight w_hom(t) in [0.0, 1.0].
        At t = 0 Mya (same species), w_hom = 1.0.
        At t >= t_crit (divergent order), w_hom drops to 0.0 to prevent false closures.
        """
        if divergence_mya <= 0.0:
            return 1.0
        if divergence_mya >= self.t_crit:
            return 0.0

        # Smooth exponential decay normalized to reach 0 at t_crit
        exp_factor = math.exp(-self.decay_lambda * divergence_mya)
        linear_gate = max(0.0, 1.0 - (divergence_mya / self.t_crit))
        return float(exp_factor * linear_gate)

    def polish_homology_with_sparse_reads(
        self,
        homology_seq: str,
        reads: List[str],
        divergence_mya: float
    ) -> Tuple[str, Dict[str, Any]]:
        """
        Uses sparse reads (e.g. 1x-3x coverage fragments) to correct evolutionary substitutions
        in the homology scaffold, restoring target-species nucleotides.
        
        Returns:
            (polished_sequence, metrics_dict)
        """
        from scripts.evaluate_metrics import needleman_wunsch

        hom = homology_seq.upper()
        if not hom:
            return "", {"substitutions_corrected": 0, "read_support_rate": 0.0}
        if not reads:
            return hom, {"substitutions_corrected": 0, "read_support_rate": 0.0}

        w_hom = self.compute_homology_weight(divergence_mya)
        n = len(hom)
        
        # Position-wise nucleotide evidence counts: {'A': count, 'C': count, ...}
        # Homology contributes w_hom prior votes to its base
        base_votes: List[Dict[str, float]] = []
        for b in hom:
            votes = {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': 0.0}
            if b in votes:
                votes[b] += w_hom
            base_votes.append(votes)

        # Align each sparse read onto the homology scaffold
        read_matches_total = 0
        read_bases_aligned = 0

        for read in reads:
            r = read.strip().upper()
            if len(r) < 3:
                continue
            
            # Global/semi-global NW alignment of read against homology scaffold
            a_hom, a_read, _ = needleman_wunsch(hom, r)
            
            h_idx = 0
            for ch, cr in zip(a_hom, a_read):
                if ch != '-':
                    if 0 <= h_idx < n:
                        if cr in ['A', 'C', 'G', 'T']:
                            # Sparse read contributes 1.0 full vote to the observed base
                            base_votes[h_idx][cr] += 1.0
                            read_bases_aligned += 1
                            if cr == ch:
                                read_matches_total += 1
                    h_idx += 1

        # Consensus selection at each position
        polished = []
        corrections = 0
        read_supported_positions = 0

        for i in range(n):
            orig_base = hom[i]
            votes = base_votes[i]
            # Find base with max votes
            best_base = max(votes.keys(), key=lambda b: votes[b])
            best_vote = votes[best_base]
            total_votes = sum(votes.values())

            # If read evidence actively votes for a different base than homology
            read_votes_for_best = votes[best_base] - (w_hom if best_base == orig_base else 0.0)
            if best_base != orig_base and read_votes_for_best >= 1.0:
                polished.append(best_base)
                corrections += 1
            else:
                polished.append(orig_base)

            if total_votes > w_hom:
                read_supported_positions += 1

        polished_seq = "".join(polished)
        support_rate = read_supported_positions / max(1, n)

        metrics = {
            "substitutions_corrected": corrections,
            "read_support_rate": float(support_rate),
            "homology_weight": float(w_hom)
        }
        return polished_seq, metrics

    def fuse(
        self,
        homology_seq: str,
        reads: List[str],
        divergence_mya: float,
        flank_identity: float,
        expected_gap_len: int
    ) -> Dict[str, Any]:
        """
        Executes full multimodal evidence fusion and assigns composite confidence score.
        """
        w_hom = self.compute_homology_weight(divergence_mya)
        polished_seq, polish_meta = self.polish_homology_with_sparse_reads(
            homology_seq=homology_seq,
            reads=reads,
            divergence_mya=divergence_mya
        )

        # Composite Confidence Scoring C(S) in [0.0, 1.0]:
        # 1. Read support factor (0.0 to 1.0)
        read_factor = polish_meta["read_support_rate"]
        # 2. Evolutionary homology factor
        hom_factor = flank_identity * w_hom
        # 3. Length discordance factor
        len_diff = abs(len(polished_seq) - expected_gap_len)
        len_factor = max(0.0, 1.0 - (len_diff / max(1, expected_gap_len)))

        # Weighted composite confidence: Read evidence (50%), Homology synteny (30%), Length (20%)
        # If no reads available, homology & length take up 60% and 40%
        if reads and len(reads) > 0:
            confidence = 0.50 * read_factor + 0.30 * hom_factor + 0.20 * len_factor
        else:
            confidence = 0.60 * hom_factor + 0.40 * len_factor

        return {
            "candidate_sequence": polished_seq,
            "confidence_score": float(confidence),
            "homology_weight": float(w_hom),
            "substitutions_corrected": polish_meta["substitutions_corrected"],
            "read_support_rate": polish_meta["read_support_rate"],
            "length_factor": float(len_factor)
        }


class SelectivePredictionGater:
    """
    Decides whether to commit a reconstructed gap sequence or reject it (retaining 'N's)
    to protect assembly integrity against False Closures.
    """
    def __init__(
        self,
        min_confidence_threshold: float = 0.70,
        max_allowed_len_error: float = 0.20,
        divergence_cutoff_mya: float = EMPIRICAL_T_CRIT
    ):
        self.min_confidence = min_confidence_threshold
        self.max_len_error = max_allowed_len_error
        self.divergence_cutoff = divergence_cutoff_mya

    def evaluate_gating(
        self,
        fusion_result: Dict[str, Any],
        divergence_mya: float,
        expected_gap_len: int,
        reads_present: bool = True
    ) -> Dict[str, Any]:
        """
        Applies strict gating rules.
        If rejected, returns 'N'*expected_gap_len.
        """
        cand = fusion_result.get("candidate_sequence", "")
        conf = fusion_result.get("confidence_score", 0.0)
        
        reasons = []

        # Rule 1: Empty or zero-length prediction
        if not cand:
            reasons.append("Empty candidate sequence")

        # Rule 2: Excessive length discordance (> 20%)
        len_error = abs(len(cand) - expected_gap_len) / max(1, expected_gap_len)
        if len_error > self.max_len_error:
            reasons.append(f"Length error exceeds {self.max_len_error*100:.0f}% ({len_error*100:.1f}%)")

        # Rule 3: Extreme evolutionary distance without read support
        if divergence_mya >= self.divergence_cutoff and not reads_present:
            reasons.append(f"Divergence {divergence_mya:.1f} Mya exceeds critical threshold {self.divergence_cutoff:.1f} Mya without read support")

        # Rule 4: Composite confidence threshold
        if conf < self.min_confidence:
            reasons.append(f"Confidence score {conf:.2f} below threshold {self.min_confidence:.2f}")

        if reasons:
            return {
                "final_sequence": "N" * expected_gap_len,
                "is_accepted": False,
                "confidence_score": conf,
                "rejection_reason": "; ".join(reasons)
            }

        return {
            "final_sequence": cand,
            "is_accepted": True,
            "confidence_score": conf,
            "rejection_reason": None
        }
