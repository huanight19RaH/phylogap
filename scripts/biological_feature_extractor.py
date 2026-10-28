"""
PhyloGap - Biological Feature Extractor & Preprocessing Engine.
Implements:
1. Data Preprocessing:
   - FlankUniquenessFilter: Evaluates k-mer frequency/entropy, expands flank dynamically if repetitive.
   - SparseReadQualityFilter: Filters low-quality/chimeric reads, extracts localized overlapping reads.
   - GapBiologicalProfiler: Analyzes local GC-isochore context, detects microsatellite/homopolymer tracts.
2. Evolutionary & Functional Feature Extraction:
   - KimuraEvolutionaryMatrix: Implements Kimura 2-Parameter (K80) transition/transversion dynamics.
   - CodonFrameExtractor: Identifies triplet reading frame and codon positions (1, 2, 3).
   - RichFeatureVectorBuilder: Constructs 22-dimensional feature tensor per nucleotide locus.

Author: PhyloGap Research Team
Date: 2026-10-28
"""

import math
from typing import List, Dict, Any, Tuple, Optional
import numpy as np

# Neutral avian mutation parameters
DEFAULT_TRANSITION_TRANSVERSION_RATIO = 2.5 # kappa = alpha / beta in birds
AVIAN_NEUTRAL_GC = 0.42

class FlankUniquenessFilter:
    """
    Evaluates flank sequence complexity and repeat density.
    Expands flank dynamically if k-mer entropy is low or repeats are detected,
    preventing ambiguous multi-mapping in genome gap closure.
    """
    def __init__(self, kmer_size: int = 5, min_entropy_threshold: float = 1.35):
        self.kmer_size = kmer_size
        self.min_entropy = min_entropy_threshold

    def calculate_kmer_entropy(self, sequence: str) -> float:
        """
        Calculates Shannon entropy of k-mers in the sequence:
        H = - sum(p * log2(p))
        Low entropy indicates simple sequence repeats (SSRs) or homopolymers.
        """
        seq = sequence.upper().strip()
        n = len(seq)
        if n < self.kmer_size:
            return 0.0

        kmers: Dict[str, int] = {}
        total_kmers = n - self.kmer_size + 1
        for i in range(total_kmers):
            k = seq[i:i + self.kmer_size]
            kmers[k] = kmers.get(k, 0) + 1

        entropy = 0.0
        for count in kmers.values():
            p = count / total_kmers
            if p > 0:
                entropy -= p * math.log2(p)
        return float(entropy)

    def detect_homopolymers(self, sequence: str, min_run_length: int = 5) -> List[Tuple[str, int, int]]:
        """
        Detects mononucleotide repeat runs (e.g. AAAAA) that trap neural models.
        Returns list of (base, start_pos, length).
        """
        seq = sequence.upper().strip()
        runs = []
        if not seq:
            return runs

        current_base = seq[0]
        current_len = 1
        start_idx = 0

        for i in range(1, len(seq)):
            if seq[i] == current_base:
                current_len += 1
            else:
                if current_len >= min_run_length:
                    runs.append((current_base, start_idx, current_len))
                current_base = seq[i]
                current_len = 1
                start_idx = i

        if current_len >= min_run_length:
            runs.append((current_base, start_idx, current_len))

        return runs

    def assess_flank_anchor(
        self,
        flank_seq: str,
        extended_context: Optional[str] = None
    ) -> Dict[str, Any]:
        """
        Audits flank anchor stability. If flank is repetitive and extended_context is supplied,
        recommends expansion length.
        """
        seq = flank_seq.upper().strip()
        entropy = self.calculate_kmer_entropy(seq)
        homopolymers = self.detect_homopolymers(seq)

        is_ambiguous = (entropy < self.min_entropy) or (len(homopolymers) > 0 and any(r[2] >= 5 for r in homopolymers))

        expanded_seq = seq
        expansion_applied = False

        if is_ambiguous and extended_context:
            ext_clean = extended_context.upper().strip()
            if len(ext_clean) > len(seq):
                expanded_seq = ext_clean
                expansion_applied = True
                entropy = self.calculate_kmer_entropy(expanded_seq)
                is_ambiguous = entropy < self.min_entropy

        return {
            "original_length": len(seq),
            "final_flank": expanded_seq,
            "kmer_entropy": float(entropy),
            "is_ambiguous": bool(is_ambiguous),
            "homopolymer_runs": homopolymers,
            "expansion_applied": bool(expansion_applied)
        }


class SparseReadQualityFilter:
    """
    Cleans, deduplicates, and bins sparse sequencing reads covering the gap.
    """
    def __init__(self, min_read_len: int = 15, max_homopolymer_fraction: float = 0.60):
        self.min_read_len = min_read_len
        self.max_homopolymer_fraction = max_homopolymer_fraction

    def filter_and_bin_reads(
        self,
        raw_reads: List[str],
        flank_left: str,
        flank_right: str,
        homology_seq: str,
        kmer_match_len: int = 12
    ) -> List[str]:
        """
        Filters out low-complexity/degenerate reads and retains only reads that share
        at least one k-mer with the gap context (flanks or homologous scaffold).
        """
        valid_bases = set(['A', 'C', 'G', 'T'])
        target_context = (flank_left[-30:] + homology_seq + flank_right[:30]).upper()
        
        # Build reference k-mer set
        ref_kmers = set()
        if len(target_context) >= kmer_match_len:
            for i in range(len(target_context) - kmer_match_len + 1):
                ref_kmers.add(target_context[i:i + kmer_match_len])

        cleaned_reads: List[str] = []
        for r in raw_reads:
            read_clean = r.upper().strip()
            if len(read_clean) < self.min_read_len:
                continue

            # Check nucleotide alphabet
            if not all(ch in valid_bases for ch in read_clean):
                continue

            # Check homopolymer saturation (e.g. AAAAAAAAA)
            max_single_base_count = max(read_clean.count(b) for b in valid_bases)
            if (max_single_base_count / len(read_clean)) > self.max_homopolymer_fraction:
                continue

            # Check localized overlap with target locus
            if ref_kmers:
                has_overlap = False
                for i in range(len(read_clean) - kmer_match_len + 1):
                    k = read_clean[i:i + kmer_match_len]
                    if k in ref_kmers:
                        has_overlap = True
                        break
                if not has_overlap:
                    continue

            cleaned_reads.append(read_clean)

        return cleaned_reads


class GapBiologicalProfiler:
    """
    Profiles local genomic environment: GC content, Isochore categorization,
    and expected evolutionary stability.
    """
    @staticmethod
    def compute_gc_fraction(sequence: str) -> float:
        seq = sequence.upper().strip()
        if not seq:
            return AVIAN_NEUTRAL_GC
        gc = sum(1 for ch in seq if ch in ('G', 'C'))
        return float(gc / len(seq))

    @staticmethod
    def categorize_isochore(gc_fraction: float) -> str:
        """
        Classifies avian isochore families:
        L1, L2 (GC < 40%): AT-rich macrochromosomes, slower mutation.
        H1, H2, H3 (GC >= 46%): GC-rich microchromosomes, high recombination/mutation.
        """
        if gc_fraction < 0.38:
            return "L1_AT_rich"
        elif gc_fraction < 0.43:
            return "L2_neutral"
        elif gc_fraction < 0.48:
            return "H1_moderate_GC"
        elif gc_fraction < 0.53:
            return "H2_elevated_GC"
        else:
            return "H3_microchromosome_hotspot"


class KimuraEvolutionaryMatrix:
    """
    Implements Kimura 2-Parameter (K80) nucleotide substitution model.
    Accounts for differing rates between Transitions (A <-> G, C <-> T)
    and Transversions (A,G <-> C,T).

    P(i -> j, t):
    - Same base (identity):
        P_ii = 1/4 + 1/4 * exp(-4*beta*t) + 1/2 * exp(-2*(alpha + beta)*t)
    - Transition (A <-> G, C <-> T):
        P_ij_ts = 1/4 + 1/4 * exp(-4*beta*t) - 1/2 * exp(-2*(alpha + beta)*t)
    - Transversion (A,G <-> C,T):
        P_ij_tv = 1/4 - 1/4 * exp(-4*beta*t)
    """
    def __init__(self, base_rate: float = 0.006395, kappa: float = DEFAULT_TRANSITION_TRANSVERSION_RATIO):
        self.base_rate = base_rate
        self.kappa = kappa # kappa = alpha / beta
        # Overall substitution rate u = (alpha + 2*beta) / 2
        # Solving for beta and alpha:
        # u = beta * (kappa + 2) / 2  =>  beta = 2 * u / (kappa + 2)
        self.beta = (2.0 * self.base_rate) / (self.kappa + 2.0)
        self.alpha = self.kappa * self.beta

    def transition_probabilities(self, divergence_mya: float) -> Dict[str, float]:
        """
        Calculates P_identical, P_transition, P_transversion at divergence time t (Mya).
        """
        t = max(0.0, float(divergence_mya))
        if t == 0.0:
            return {"p_ident": 1.0, "p_ts": 0.0, "p_tv": 0.0}

        exp_4beta = math.exp(-4.0 * self.beta * t)
        exp_2alphabeta = math.exp(-2.0 * (self.alpha + self.beta) * t)

        p_ident = 0.25 + 0.25 * exp_4beta + 0.50 * exp_2alphabeta
        p_ts = 0.25 + 0.25 * exp_4beta - 0.50 * exp_2alphabeta
        p_tv = 0.25 - 0.25 * exp_4beta

        # Normalize to ensure sum equals 1.0 (1 ident + 1 ts + 2 tv = 1.0)
        total = p_ident + p_ts + 2.0 * p_tv
        return {
            "p_ident": max(0.0, float(p_ident / total)),
            "p_ts": max(0.0, float(p_ts / total)),
            "p_tv": max(0.0, float(p_tv / total))
        }

    def get_substitution_probability(self, source_base: str, target_base: str, divergence_mya: float) -> float:
        """
        Returns conditional probability P(target_base | source_base, t).
        """
        src = source_base.upper()
        tgt = target_base.upper()
        probs = self.transition_probabilities(divergence_mya)

        if src == tgt:
            return probs["p_ident"]

        # Transitions: A <-> G (purines) or C <-> T (pyrimidines)
        is_transition = (
            (src == 'A' and tgt == 'G') or (src == 'G' and tgt == 'A') or
            (src == 'C' and tgt == 'T') or (src == 'T' and tgt == 'C')
        )
        if is_transition:
            return probs["p_ts"]
        else:
            return probs["p_tv"]


class CodonFrameExtractor:
    """
    Handles triplet codon structure in coding exon regions.
    Tracks reading frame index (1, 2, 3) to enforce purifying selection on pos 1 & 2.
    """
    @staticmethod
    def get_codon_position(locus_idx: int, frame_offset: int = 0) -> int:
        """
        Returns codon position: 1, 2, or 3 (wobble).
        """
        return ((locus_idx + frame_offset) % 3) + 1

    @staticmethod
    def is_reading_frame_preserved(gap_length: int) -> bool:
        """
        In coding exons, indels must be multiples of 3 bp to avoid catastrophic frameshifts.
        """
        return (gap_length % 3) == 0


class RichFeatureVectorBuilder:
    """
    Constructs 22-dimensional rich biological feature tensor per nucleotide locus.
    """
    def __init__(self, k80_model: Optional[KimuraEvolutionaryMatrix] = None):
        self.k80 = k80_model or KimuraEvolutionaryMatrix()

    def build_locus_feature_vector(
        self,
        homology_base: str,
        read_counts: Dict[str, float],
        divergence_mya: float,
        feature_type: str = "intron",
        locus_idx: int = 0,
        frame_offset: int = 0,
        local_gc_fraction: float = AVIAN_NEUTRAL_GC,
        homology_weight: float = 1.0
    ) -> np.ndarray:
        """
        Vector schema (22 dimensions):
        - [0:4]   One-hot Homology base: [A, C, G, T]
        - [4:8]   Kimura K80 mutation prior: [P(A|Hom,t), P(C|Hom,t), P(G|Hom,t), P(T|Hom,t)]
        - [8:12]  Sparse reads count normalized: [c_A, c_C, c_G, c_T]
        - [12:15] One-hot Genomic feature: [Exon, Intron, Repeat]
        - [15:19] One-hot Codon position: [None/Non-coding, Pos_1, Pos_2, Pos_3]
        - [19]    Local GC fraction (float in [0, 1])
        - [20]    Divergence time t normalized (t / 100.0)
        - [21]    Homology weight w_hom (float in [0, 1])
        """
        vec = np.zeros(22, dtype=np.float32)
        bases = ['A', 'C', 'G', 'T']
        hom = homology_base.upper()

        # 1. Homology One-hot
        if hom in bases:
            vec[bases.index(hom)] = 1.0

        # 2. Kimura K80 Mutation Prior
        if hom in bases:
            for i, b in enumerate(bases):
                vec[4 + i] = self.k80.get_substitution_probability(hom, b, divergence_mya)

        # 3. Read counts
        total_reads = sum(read_counts.get(b, 0.0) for b in bases)
        denom = max(1.0, total_reads)
        for i, b in enumerate(bases):
            vec[8 + i] = float(read_counts.get(b, 0.0) / denom)

        # 4. Genomic feature type
        feat = feature_type.lower()
        if "exon" in feat or "coding" in feat:
            vec[12] = 1.0
        elif "repeat" in feat:
            vec[14] = 1.0
        else: # intron / intergenic
            vec[13] = 1.0

        # 5. Codon position
        if "exon" in feat or "coding" in feat:
            pos = CodonFrameExtractor.get_codon_position(locus_idx, frame_offset)
            vec[15 + pos] = 1.0
        else:
            vec[15] = 1.0 # None / Non-coding

        # 6. Isochore GC
        vec[19] = float(local_gc_fraction)

        # 7. Divergence time t
        vec[20] = float(min(1.0, divergence_mya / 100.0))

        # 8. Homology weight
        vec[21] = float(homology_weight)

        return vec
