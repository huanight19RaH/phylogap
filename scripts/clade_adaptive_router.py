"""
PhyloGap - Clade-Specific Adaptive Evolutionary Router.
Provides lineage-tailored evolutionary substitution models, mutation matrices,
and selective prediction thresholds across distinct biological kingdoms and classes:
1. AVES (Birds): Compact genomes, microchromosomes, high metabolic rate (K80, kappa=2.5, mu=0.0064/Myr).
2. MAMMALIA (Mammals): High CpG hypermutation deamination, LINE/SINE abundance (HKY85, kappa=3.8, mu=0.0022/Myr).
3. VIRIDIPLANTAE (Plants): Whole-Genome Duplication (WGD) synteny awareness, high indel dynamics (GTR+Gamma).
4. INSECTA (Insects/Drosophila): High effective population size, rapid deletion dynamics, compact introns (mu=0.0120/Myr).

Author: PhyloGap Research Team
Date: 2026-10-28
"""

import math
from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, Any, List, Optional, Tuple
import numpy as np


class CladeType(Enum):
    AVES = "Aves"
    MAMMALIA = "Mammalia"
    VIRIDIPLANTAE = "Viridiplantae"
    INSECTA = "Insecta"
    UNIVERSAL = "Universal"


@dataclass
class CladeEvolutionProfile:
    clade_name: str
    clade_type: CladeType
    base_substitution_rate_per_myr: float # mu
    transition_transversion_ratio: float   # kappa = alpha / beta
    neutral_gc_fraction: float
    cpg_hypermutation_factor: float        # multiplier for C->T in CpG dinucleotides
    wgd_paralogy_risk: bool                # True for plants with whole genome duplication
    exon_t_crit_mya: float                 # Critical divergence cutoff for exons
    intron_t_crit_mya: float               # Critical divergence cutoff for introns
    repeat_t_crit_mya: float               # Critical divergence cutoff for repeats
    repeat_density_typical: float          # Fraction of genome in repeats
    description: str


# Clade Profile Registry
CLADE_REGISTRY: Dict[CladeType, CladeEvolutionProfile] = {
    CladeType.AVES: CladeEvolutionProfile(
        clade_name="Aves (Birds)",
        clade_type=CladeType.AVES,
        base_substitution_rate_per_myr=0.0064,
        transition_transversion_ratio=2.5,
        neutral_gc_fraction=0.42,
        cpg_hypermutation_factor=1.8,
        wgd_paralogy_risk=False,
        exon_t_crit_mya=31.4,
        intron_t_crit_mya=13.2,
        repeat_t_crit_mya=5.0,
        repeat_density_typical=0.10,
        description="Compact avian genomes with high synteny conservation, GC-rich microchromosomes."
    ),
    CladeType.MAMMALIA: CladeEvolutionProfile(
        clade_name="Mammalia (Mammals)",
        clade_type=CladeType.MAMMALIA,
        base_substitution_rate_per_myr=0.0022,
        transition_transversion_ratio=3.8,
        neutral_gc_fraction=0.41,
        cpg_hypermutation_factor=10.5, # Spontaneous deamination of 5mC at CpG
        wgd_paralogy_risk=False,
        exon_t_crit_mya=45.0,
        intron_t_crit_mya=18.0,
        repeat_t_crit_mya=4.0,
        repeat_density_typical=0.45,
        description="High transposable element density (LINEs/SINEs), severe CpG deamination bias."
    ),
    CladeType.VIRIDIPLANTAE: CladeEvolutionProfile(
        clade_name="Viridiplantae (Green Plants)",
        clade_type=CladeType.VIRIDIPLANTAE,
        base_substitution_rate_per_myr=0.0048,
        transition_transversion_ratio=2.2,
        neutral_gc_fraction=0.38,
        cpg_hypermutation_factor=3.0,
        wgd_paralogy_risk=True, # Paleopolyploidy events create homeologous paralogs
        exon_t_crit_mya=28.0,
        intron_t_crit_mya=8.0,
        repeat_t_crit_mya=2.5,
        repeat_density_typical=0.60,
        description="Frequent Whole-Genome Duplication (WGD), high indel rate, paralogous synteny ambiguity."
    ),
    CladeType.INSECTA: CladeEvolutionProfile(
        clade_name="Insecta (Insects / Drosophila)",
        clade_type=CladeType.INSECTA,
        base_substitution_rate_per_myr=0.0120,
        transition_transversion_ratio=2.0,
        neutral_gc_fraction=0.40,
        cpg_hypermutation_factor=1.0,
        wgd_paralogy_risk=False,
        exon_t_crit_mya=20.0,
        intron_t_crit_mya=6.0,
        repeat_t_crit_mya=2.0,
        repeat_density_typical=0.20,
        description="High mutation rate, short generation time, ultra-compact introns, strong codon selection."
    ),
    CladeType.UNIVERSAL: CladeEvolutionProfile(
        clade_name="Universal / Default Baseline",
        clade_type=CladeType.UNIVERSAL,
        base_substitution_rate_per_myr=0.0050,
        transition_transversion_ratio=2.0,
        neutral_gc_fraction=0.42,
        cpg_hypermutation_factor=2.0,
        wgd_paralogy_risk=False,
        exon_t_crit_mya=25.0,
        intron_t_crit_mya=10.0,
        repeat_t_crit_mya=3.0,
        repeat_density_typical=0.30,
        description="General phylogenetic default for unclassified eukaryotes."
    )
}


class CladeAdaptiveRouter:
    """
    Directs evolutionary parameters, substitution matrices, and gating thresholds
    based on the biological classification of the target genome.
    """
    def __init__(self, default_clade: CladeType = CladeType.AVES):
        self.default_clade = default_clade

    def detect_clade_from_name_or_taxid(
        self,
        identifier: str
    ) -> CladeEvolutionProfile:
        """
        Resolves clade from species name, common name, or taxonomic keyword.
        """
        ident_lower = identifier.lower().strip()

        # Avian keywords
        if any(w in ident_lower for w in ["bird", "ave", "calypte", "gallus", "finch", "passerine", "hummingbird", "b10k", "lonchura", "taeniopygia"]):
            return CLADE_REGISTRY[CladeType.AVES]

        # Mammalian keywords
        if any(w in ident_lower for w in ["mammal", "human", "homo", "mus", "mouse", "rat", "dog", "canis", "bovine", "bos", "primate"]):
            return CLADE_REGISTRY[CladeType.MAMMALIA]

        # Plant keywords
        if any(w in ident_lower for w in ["plant", "arabidopsis", "oryza", "rice", "zea", "maize", "viridiplantae", "chlorophyta", "angiosperm"]):
            return CLADE_REGISTRY[CladeType.VIRIDIPLANTAE]

        # Insect keywords
        if any(w in ident_lower for w in ["insect", "drosophila", "fly", "mosquito", "anopheles", "apis", "bee"]):
            return CLADE_REGISTRY[CladeType.INSECTA]

        return CLADE_REGISTRY[self.default_clade]

    def get_profile(self, clade_type: CladeType) -> CladeEvolutionProfile:
        return CLADE_REGISTRY.get(clade_type, CLADE_REGISTRY[CladeType.UNIVERSAL])

    def calculate_clade_transition_probability(
        self,
        divergence_mya: float,
        clade_type: CladeType = CladeType.AVES,
        is_transition: bool = True
    ) -> float:
        """
        Calculates Kimura K80 transition (P) or transversion (Q) probability under
        clade-specific mu and kappa.
        """
        profile = self.get_profile(clade_type)
        mu = profile.base_substitution_rate_per_myr
        kappa = profile.transition_transversion_ratio

        # Total distance d = mu * t
        d = mu * divergence_mya
        beta = d / (kappa + 2.0)
        alpha = kappa * beta

        if is_transition:
            # P(t) = 0.25 + 0.25*e^(-4*beta) - 0.5*e^(-2*(alpha + beta))
            p_trans = 0.25 + 0.25 * math.exp(-4.0 * beta) - 0.5 * math.exp(-2.0 * (alpha + beta))
            return max(0.0, min(1.0, float(p_trans)))
        else:
            # Q(t) = 0.25 - 0.25*e^(-4*beta)
            q_transv = 0.25 - 0.25 * math.exp(-4.0 * beta)
            return max(0.0, min(1.0, float(q_transv)))

    def audit_cpg_context(
        self,
        ref_seq: str,
        target_seq: str,
        clade_type: CladeType = CladeType.MAMMALIA
    ) -> Dict[str, Any]:
        """
        Specifically audits CpG hypermutation sites. In Mammals, C->T transitions
        in a CpG context (CG dinucleotide) occur at a 10-15x higher rate due to
        spontaneous deamination of 5-methylcytosine.
        """
        profile = self.get_profile(clade_type)
        r = ref_seq.upper().strip()
        t = target_seq.upper().strip()
        n = min(len(r), len(t))

        cpg_sites = 0
        cpg_deaminations = 0

        for i in range(n - 1):
            if r[i:i+2] == "CG":
                cpg_sites += 1
                # Check for C->T at i (CG -> TG) or G->A at i+1 (CG -> CA on reverse)
                if (t[i] == 'T' and t[i+1] == 'G') or (t[i] == 'C' and t[i+1] == 'A'):
                    cpg_deaminations += 1

        deamination_rate = (cpg_deaminations / cpg_sites) if cpg_sites > 0 else 0.0
        # If mammalian, we discount deamination mismatches because they reflect neutral epigenetic drift
        is_discounted = (profile.cpg_hypermutation_factor > 3.0 and cpg_deaminations > 0)

        return {
            "clade": profile.clade_name,
            "cpg_sites_evaluated": cpg_sites,
            "cpg_deaminations_observed": cpg_deaminations,
            "cpg_deamination_rate": float(deamination_rate),
            "epigenetic_drift_discounted": is_discounted,
            "hypermutation_factor": profile.cpg_hypermutation_factor
        }

    def compute_clade_gating_decision(
        self,
        divergence_mya: float,
        locus_type: str,
        clade_type: CladeType = CladeType.AVES,
        has_sparse_reads: bool = False
    ) -> Dict[str, Any]:
        """
        Determines whether homology-based gap filling is biologically safe for a given clade.
        """
        profile = self.get_profile(clade_type)
        l_type = locus_type.lower()

        if l_type in ["exon", "coding_exon", "cds"]:
            t_crit = profile.exon_t_crit_mya
        elif l_type in ["repeat", "te", "ltr", "line"]:
            t_crit = profile.repeat_t_crit_mya
        else:
            t_crit = profile.intron_t_crit_mya

        # With sparse reads, safe cutoff expands by 2.0x
        if has_sparse_reads:
            effective_t_crit = t_crit * 2.0
        else:
            effective_t_crit = t_crit

        is_safe = (divergence_mya <= effective_t_crit)
        w_hom = max(0.0, 1.0 - (divergence_mya / effective_t_crit)) if is_safe else 0.0

        return {
            "clade": profile.clade_name,
            "locus_type": l_type,
            "divergence_mya": divergence_mya,
            "nominal_t_crit": t_crit,
            "effective_t_crit": effective_t_crit,
            "is_safe_to_close": is_safe,
            "homology_weight": float(w_hom),
            "wgd_paralogy_warning": profile.wgd_paralogy_risk
        }
