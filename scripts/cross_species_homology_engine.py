"""
PhyloGap - Cross-Species Homology Gap Reconstruction Engine.
Implementation of Synteny-Aware Homologous Anchor Mapping,
Candidate Sequence Extraction, and Selective Prediction Gating.

Addresses Research Question RQ3:
"How does evolutionary divergence between target and reference species
correlate with nucleotide-level reconstruction accuracy and False-Closure Rate?"
"""

import math
import random
from typing import Dict, Any, Tuple, Optional, List
import numpy as np

# Global random seed for reproducible benchmark simulations
RANDOM_SEED = 42

class SyntenicAnchorMapper:
    """
    Identifies syntenic anchors between flanking sequences of a target gap
    and a homologous reference genome locus.
    """
    def __init__(self, kmer_size: int = 15, max_flank_mismatch: float = 0.25):
        self.kmer_size = kmer_size
        self.max_flank_mismatch = max_flank_mismatch

    @staticmethod
    def reverse_complement(seq: str) -> str:
        comp = str.maketrans("ACGTUacgtuNn-", "TGCAAtgcaaNn-")
        return seq.translate(comp)[::-1]

    def _find_flank_pos(self, query: str, target: str, k_sizes: List[int] = [11, 8, 6]) -> Tuple[int, float]:
        """
        Locates query inside target using multi-scale k-mer voting and bounded Needleman-Wunsch.
        Returns: (best_pos, identity_ratio)
        """
        from scripts.evaluate_metrics import needleman_wunsch
        q = query.upper()
        t = target.upper()
        len_q, len_t = len(q), len(t)
        if len_q == 0 or len_t == 0:
            return -1, 0.0

        best_pos = -1
        for k in k_sizes:
            if len_q < k:
                continue
            votes = {}
            step = max(1, k // 3)
            for i in range(0, len_q - k + 1, step):
                kmer = q[i:i + k]
                idx = 0
                while True:
                    idx = t.find(kmer, idx)
                    if idx == -1:
                        break
                    pos = idx - i
                    if 0 <= pos <= len_t - len_q + 50:
                        votes[pos] = votes.get(pos, 0) + 1
                    idx += 1
            if votes:
                best_pos = max(votes.keys(), key=lambda p: votes[p])
                break

        if best_pos == -1:
            # Fallback: sliding window
            best_score = -1
            for p in range(0, max(1, len_t - len_q + 1), max(1, len_q // 10)):
                sub = t[p:p + len_q]
                m = sum(1 for c1, c2 in zip(q, sub) if c1 == c2 and c1 != 'N')
                if m > best_score:
                    best_score = m
                    best_pos = p

        if best_pos < 0 or best_pos >= len_t:
            return -1, 0.0

        # Identity computation with Needleman-Wunsch against target window
        sub_len = min(len_t - best_pos, len_q)
        sub = t[best_pos:best_pos + sub_len]
        a_q, a_t, _ = needleman_wunsch(q, sub)
        matches = sum(1 for cq, ct in zip(a_q, a_t) if cq == ct and cq != '-' and cq != 'N')
        id_ratio = matches / max(1, len_q)
        return best_pos, float(id_ratio)

    def map_flanks_to_reference(
        self,
        flank_left: str,
        flank_right: str,
        reference_locus: str,
        expected_gap_len: int,
        search_radius_factor: float = 3.0
    ) -> Dict[str, Any]:
        """
        Maps left and right flanks to reference locus, verifying syntenic collinearity.
        """
        flank_l = flank_left.upper()
        flank_r = flank_right.upper()
        ref = reference_locus.upper()

        if len(flank_l) == 0 or len(flank_r) == 0 or len(ref) == 0:
            return {
                "is_mapped": False,
                "orientation": "unmapped",
                "extracted_gap": "",
                "left_identity": 0.0,
                "right_identity": 0.0,
                "length_discordance": 1.0,
                "error": "Empty input sequence"
            }

        # 1. Forward alignment test
        pos_l, id_l = self._find_flank_pos(flank_l, ref)
        pos_r, id_r = self._find_flank_pos(flank_r, ref)

        if pos_l != -1 and pos_r != -1:
            left_end = pos_l + len(flank_l)
            right_start = pos_r

            # Verify synteny: right flank must come after left flank
            if right_start > left_end:
                extracted = ref[left_end:right_start]
                len_diff = abs(len(extracted) - expected_gap_len)
                rel_len_diff = len_diff / max(1, expected_gap_len)

                # Check if flanks pass minimum identity
                if id_l >= (1.0 - self.max_flank_mismatch) and id_r >= (1.0 - self.max_flank_mismatch):
                    return {
                        "is_mapped": True,
                        "orientation": "forward",
                        "left_end_pos": left_end,
                        "right_start_pos": right_start,
                        "extracted_gap": extracted,
                        "left_identity": id_l,
                        "right_identity": id_r,
                        "length_discordance": rel_len_diff
                    }

        # 2. Check reverse complement orientation (Inversion / Strand difference)
        flank_l_rc = self.reverse_complement(flank_r)
        flank_r_rc = self.reverse_complement(flank_l)

        pos_l_rc, id_l_rc = self._find_flank_pos(flank_l_rc, ref)
        pos_r_rc, id_r_rc = self._find_flank_pos(flank_r_rc, ref)

        if pos_l_rc != -1 and pos_r_rc != -1:
            left_end_rc = pos_l_rc + len(flank_l_rc)
            right_start_rc = pos_r_rc

            if right_start_rc > left_end_rc:
                extracted_rc = ref[left_end_rc:right_start_rc]
                extracted_fwd = self.reverse_complement(extracted_rc)
                len_diff = abs(len(extracted_fwd) - expected_gap_len)
                rel_len_diff = len_diff / max(1, expected_gap_len)

                if id_l_rc >= (1.0 - self.max_flank_mismatch) and id_r_rc >= (1.0 - self.max_flank_mismatch):
                    return {
                        "is_mapped": True,
                        "orientation": "reverse_complement",
                        "left_end_pos": left_end_rc,
                        "right_start_pos": right_start_rc,
                        "extracted_gap": extracted_fwd,
                        "left_identity": id_l_rc,
                        "right_identity": id_r_rc,
                        "length_discordance": rel_len_diff
                    }

        return {
            "is_mapped": False,
            "orientation": "discordant" if (id_l > 0.5 or id_r > 0.5) else "unmapped",
            "extracted_gap": "",
            "left_identity": max(id_l, id_l_rc if 'id_l_rc' in locals() else 0.0),
            "right_identity": max(id_r, id_r_rc if 'id_r_rc' in locals() else 0.0),
            "length_discordance": 1.0,
            "error": "Flanks could not be confidently anchored within syntenic bounds"
        }

        return {
            "is_mapped": False,
            "orientation": "discordant" if (id_l > 0.5 or id_r > 0.5) else "unmapped",
            "extracted_gap": "",
            "left_identity": max(id_l, id_l_rc),
            "right_identity": max(id_r, id_r_rc),
            "length_discordance": 1.0,
            "error": "Flanks could not be confidently anchored within syntenic bounds"
        }


class EvidenceAwareSelectivePredictor:
    """
    Evaluates candidate homology sequences and selectively accepts or rejects predictions
    to protect genome assemblies against False Closures.
    """
    def __init__(
        self,
        min_flank_identity: float = 0.82,
        max_allowed_len_discordance: float = 0.35,
        min_confidence_threshold: float = 0.70
    ):
        self.min_flank_identity = min_flank_identity
        self.max_allowed_len_discordance = max_allowed_len_discordance
        self.min_confidence_threshold = min_confidence_threshold

    def evaluate_and_filter(
        self,
        mapping_result: Dict[str, Any],
        expected_gap_len: int
    ) -> Dict[str, Any]:
        """
        Applies gating criteria to candidate gap sequence.
        
        Returns:
            Dict containing:
            - final_prediction: str (reconstructed nucleotide sequence OR 'N'*gap_len)
            - is_accepted: bool (True: accepted, False: rejected/preserved as N)
            - confidence_score: float [0.0, 1.0]
            - rejection_reason: str or None
        """
        if not mapping_result.get("is_mapped", False):
            return {
                "final_prediction": "N" * expected_gap_len,
                "is_accepted": False,
                "confidence_score": 0.0,
                "rejection_reason": "Flank synteny unmapped"
            }

        id_l = mapping_result["left_identity"]
        id_r = mapping_result["right_identity"]
        len_disc = mapping_result["length_discordance"]
        extracted = mapping_result["extracted_gap"]

        # Compute composite confidence score
        # Flank homology contributes 60%, length agreement contributes 40%
        flank_score = (id_l + id_r) / 2.0
        len_score = max(0.0, 1.0 - (len_disc / self.max_allowed_len_discordance))
        confidence = 0.60 * flank_score + 0.40 * len_score

        # Check rejection rules
        rejection_reasons = []
        if flank_score < self.min_flank_identity:
            rejection_reasons.append(f"Low flank identity ({flank_score:.2f} < {self.min_flank_identity})")
        if len_disc > self.max_allowed_len_discordance:
            rejection_reasons.append(f"Excessive length discordance ({len_disc:.2f} > {self.max_allowed_len_discordance})")
        if confidence < self.min_confidence_threshold:
            rejection_reasons.append(f"Composite confidence below threshold ({confidence:.2f} < {self.min_confidence_threshold})")

        if rejection_reasons:
            return {
                "final_prediction": "N" * expected_gap_len,
                "is_accepted": False,
                "confidence_score": confidence,
                "rejection_reason": "; ".join(rejection_reasons)
            }

        return {
            "final_prediction": extracted,
            "is_accepted": True,
            "confidence_score": confidence,
            "rejection_reason": None
        }


class AvianEvolutionarySimulator:
    """
    Biologically calibrated evolutionary divergence simulator for avian genomes.
    Implements Kimura 2-parameter (K2P) nucleotide substitution and affine indels
    based on avian macro-evolutionary rates (B10K, Jarvis et al., Stiller et al. 2024).
    """
    def __init__(
        self,
        substitution_rate_per_myr: float = 0.0019, # ~1.9 x 10^-9 sub/site/year
        ts_tv_ratio: float = 2.5,                  # Transition/Transversion ratio in birds
        indel_rate_ratio: float = 0.15,            # Indel events relative to substitutions
        seed: int = RANDOM_SEED
    ):
        self.sub_rate = substitution_rate_per_myr
        self.ts_tv_ratio = ts_tv_ratio
        self.indel_rate = substitution_rate_per_myr * indel_rate_ratio
        self.rng = random.Random(seed)
        self.np_rng = np.random.RandomState(seed)

    def evolve_sequence(self, sequence: str, divergence_time_mya: float) -> Tuple[str, Dict[str, int]]:
        """
        Simulates neutral evolution of a sequence over divergence_time_mya (million years).
        
        Returns:
            (evolved_sequence, mutation_stats)
        """
        seq = sequence.upper()
        if not seq:
            return "", {"substitutions": 0, "insertions": 0, "deletions": 0}

        # K2P substitution probability: p = 1 - exp(-2 * mu * t)
        # Note: 2*t because both lineages diverge from common ancestor
        total_time = 2.0 * divergence_time_mya
        p_sub = 1.0 - math.exp(-self.sub_rate * total_time)
        p_indel = 1.0 - math.exp(-self.indel_rate * total_time)

        transitions = {'A': 'G', 'G': 'A', 'C': 'T', 'T': 'C'}
        transversions = {
            'A': ['C', 'T'],
            'G': ['C', 'T'],
            'C': ['A', 'G'],
            'T': ['A', 'G']
        }
        bases = ['A', 'C', 'G', 'T']

        evolved = []
        n_sub, n_ins, n_del = 0, 0, 0

        i = 0
        while i < len(seq):
            base = seq[i]
            if base not in transitions:
                evolved.append(base)
                i += 1
                continue

            # 1. Indel check
            if self.rng.random() < p_indel:
                if self.rng.random() < 0.5:
                    # Deletion: geometric length (mean 2 bp)
                    del_len = min(len(seq) - i, self.np_rng.geometric(0.5))
                    n_del += del_len
                    i += del_len
                    continue
                else:
                    # Insertion: geometric length
                    ins_len = self.np_rng.geometric(0.5)
                    inserted_bases = "".join(self.rng.choices(bases, k=ins_len))
                    evolved.append(inserted_bases)
                    n_ins += ins_len

            # 2. Substitution check
            if self.rng.random() < p_sub:
                n_sub += 1
                # Transition vs Transversion
                prob_ts = self.ts_tv_ratio / (self.ts_tv_ratio + 1.0)
                if self.rng.random() < prob_ts:
                    evolved.append(transitions[base])
                else:
                    evolved.append(self.rng.choice(transversions[base]))
            else:
                evolved.append(base)

            i += 1

        evolved_seq = "".join(evolved)
        stats = {
            "substitutions": n_sub,
            "insertions": n_ins,
            "deletions": n_del
        }
        return evolved_seq, stats
