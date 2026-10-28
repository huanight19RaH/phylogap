"""
PhyloGap Interactive Web Dashboard - Backend Server.
FastAPI backend providing real-time evidence-aware gap reconstruction,
evolutionary risk calculation, and benchmark exploration.
"""

import os
import sys
import math
from typing import List, Optional, Dict, Any
from fastapi import FastAPI, HTTPException
from fastapi.staticfiles import StaticFiles
from fastapi.responses import FileResponse
from pydantic import BaseModel

# Ensure project root is in sys.path
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from scripts.phylogap_cli import reconstruct_gap_pipeline, DEMO_PRESETS
from scripts.cross_species_homology_engine import AvianEvolutionarySimulator
from scripts.evidence_fusion_engine import (
    EMPIRICAL_DECAY_LAMBDA,
    EMPIRICAL_T_CRIT,
    EvidenceWeightedFuser
)

app = FastAPI(
    title="PhyloGap Studio",
    description="Evidence-Aware Genome Assembly Gap Reconstruction API",
    version="1.0.0"
)

STATIC_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "static")

class ReconstructionRequest(BaseModel):
    flank_left: str
    flank_right: str
    homology_reference: str
    divergence_mya: float
    expected_gap_len: int
    sparse_reads: Optional[List[str]] = []
    ground_truth_gap: Optional[str] = None
    feature_type: Optional[str] = "intron"
    gc_fraction: Optional[float] = None
    enable_preprocessing: Optional[bool] = True

@app.get("/api/health")
def health_check():
    return {
        "status": "healthy",
        "project": "PhyloGap",
        "version": "1.0.0",
        "decay_lambda": EMPIRICAL_DECAY_LAMBDA,
        "t_crit_mya": EMPIRICAL_T_CRIT
    }

@app.get("/api/presets")
def get_presets():
    """Returns preloaded validated demo cases for 1-click loading."""
    presets = []
    for tier_id, data in DEMO_PRESETS.items():
        # Generate simulated reference if needed
        sim = AvianEvolutionarySimulator(seed=2026 + tier_id)
        target = data["flank_left"] + data["gap_truth"] + data["flank_right"]
        evolved, _ = sim.evolve_sequence(target, divergence_time_mya=data["divergence_mya"])
        full_ref = ("NNNN" * 5) + evolved + ("NNNN" * 5)

        presets.append({
            "tier_id": tier_id,
            "title": data["tier"],
            "target_species": data["target_species"],
            "reference_species": data["reference_species"],
            "divergence_mya": data["divergence_mya"],
            "flank_left": data["flank_left"],
            "gap_truth": data["gap_truth"],
            "flank_right": data["flank_right"],
            "homology_reference": full_ref,
            "sparse_reads": data["sparse_reads"],
            "expected_gap_len": len(data["gap_truth"])
        })
    return {"presets": presets}

@app.post("/api/reconstruct")
def reconstruct_gap(req: ReconstructionRequest):
    """Executes full multimodal evidence fusion & selective gating pipeline."""
    flank_l = req.flank_left.strip().upper()
    flank_r = req.flank_right.strip().upper()
    ref = req.homology_reference.strip().upper()

    if len(flank_l) < 15 or len(flank_r) < 15:
        raise HTTPException(status_code=400, detail="Flanks must be at least 15 bp for synteny anchoring.")

    if not ref:
        raise HTTPException(status_code=400, detail="Homologous reference sequence is required.")

    res = reconstruct_gap_pipeline(
        flank_left=flank_l,
        flank_right=flank_r,
        homology_reference=ref,
        divergence_mya=req.divergence_mya,
        expected_gap_len=req.expected_gap_len,
        sparse_reads=req.sparse_reads,
        ground_truth_gap=req.ground_truth_gap,
        feature_type=req.feature_type or "intron",
        gc_fraction=req.gc_fraction,
        enable_preprocessing=req.enable_preprocessing if req.enable_preprocessing is not None else True
    )

    # Base-by-base provenance decomposition for interactive visual rendering
    final_seq = res["final_sequence"]
    provenance = []
    
    if res["is_accepted"]:
        raw_hom = res.get("raw_homology_gap", "")
        for idx, base in enumerate(final_seq):
            if idx < len(raw_hom):
                orig_hom_base = raw_hom[idx]
                if base != orig_hom_base:
                    provenance.append({"base": base, "source": "read_corrected", "original": orig_hom_base})
                else:
                    provenance.append({"base": base, "source": "homology", "original": orig_hom_base})
            else:
                provenance.append({"base": base, "source": "extended", "original": "-"})
    else:
        for base in final_seq:
            provenance.append({"base": base, "source": "rejected_masked", "original": "N"})

    res["provenance"] = provenance
    return res

@app.get("/api/risk-calc")
def calculate_evolutionary_risk(mya: float, lang: str = "en"):
    """Calculates evolutionary sequence identity decay and risk tier for given Mya in English or Vietnamese."""
    t = max(0.0, float(mya))
    expected_identity = math.exp(-EMPIRICAL_DECAY_LAMBDA * t) * 100.0

    fuser = EvidenceWeightedFuser()
    w_hom = fuser.compute_homology_weight(t)

    is_vi = (lang.lower() == "vi")

    if t <= 6.0:
        risk_tier = "SAFE_CONGENERIC"
        risk_label = "Vùng An Toàn (Cùng chi / Sister Species)" if is_vi else "Safe Zone (Congeneric / Sister Species)"
        risk_color = "#10B981" # Green
        guidance = (
            "Độ chính xác tương đồng >95%. Có thể tái tạo an toàn ngay cả khi không có reads hỗ trợ."
            if is_vi else
            "Homology reference provides >95% accuracy. Safe for direct gap reconstruction even without read support."
        )
    elif t <= EMPIRICAL_T_CRIT:
        risk_tier = "CAUTION_FAMILY"
        risk_label = "Vùng Cảnh Báo (Cùng họ / Family Level)" if is_vi else "Caution Zone (Family Level Clade)"
        risk_color = "#F59E0B" # Yellow
        guidance = (
            "Độ chính xác tương đồng đạt 92-95% nhưng xuất hiện ~15% False Closure. Bắt buộc có Sparse Reads để sửa các đột biến đơn loài."
            if is_vi else
            "Homology identity is 92-95% but exhibits ~15% False Closure rate. Sparse reads required to polish species-specific substitutions."
        )
    else:
        risk_tier = "DANGER_ORDER_CLASS"
        risk_label = "Vùng Nguy Hiểm (Khác bộ / Distant Orders)" if is_vi else "Danger Zone (Distant Orders / Classes)"
        risk_color = "#F43F5E" # Red
        guidance = (
            f"Khoảng cách phân kỳ ({t:.1f} Mya) vượt mốc tới hạn t_crit ({EMPIRICAL_T_CRIT} Mya). Tỷ lệ đóng sai >95%! Bắt buộc kích hoạt Selective Gating để từ chối."
            if is_vi else
            f"Phylogenetic distance ({t:.1f} Mya) exceeds critical threshold t_crit ({EMPIRICAL_T_CRIT} Mya). False-Closure Rate >95%! Selective Gating required to reject candidate."
        )

    return {
        "divergence_mya": t,
        "expected_nucleotide_identity": round(expected_identity, 2),
        "homology_weight_w_hom": round(w_hom, 4),
        "critical_threshold_t_crit": EMPIRICAL_T_CRIT,
        "is_exceeding_t_crit": bool(t >= EMPIRICAL_T_CRIT),
        "risk_tier": risk_tier,
        "risk_label": risk_label,
        "risk_color": risk_color,
        "guidance": guidance,
        "lang": lang
    }

GENETIC_CODE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
}

SPECIES_CATALOG = [
    {
        "id": "c_anna",
        "name": "Calypte anna",
        "common_name": "Anna's Hummingbird",
        "clade": "Aves (Apodiformes)",
        "divergence_mya": 0.0,
        "genome_size_gb": 1.06,
        "gc_percent": 41.5,
        "morphology_type": "hummingbird",
        "sample_sequence": "ATGGGTCGCAAGCTGGACCTGCGCGAGTACGACCAGATCCTGGAGCGCATCGAGGAGATGCTGCAGCGCAAGGAGCTG",
        "description": "High metabolic rate, wing hovering kinematics, rapid nucleotide substitution."
    },
    {
        "id": "t_guttata",
        "name": "Taeniopygia guttata",
        "common_name": "Zebra Finch",
        "clade": "Aves (Passeriformes / Estrildidae)",
        "divergence_mya": 11.5,
        "genome_size_gb": 1.05,
        "gc_percent": 42.0,
        "morphology_type": "finch",
        "sample_sequence": "ATGGGTCGCAAGCTGGACCTGCGCGAGTATGACCAGATTCTGGAGCGCATCGAGGAGATGCTCCAGCGCAAGGAGCTG",
        "description": "Avian model organism for vocal learning, compact genome, high synteny."
    },
    {
        "id": "g_gallus",
        "name": "Gallus gallus",
        "common_name": "Red Junglefowl / Chicken",
        "clade": "Aves (Galliformes)",
        "divergence_mya": 88.0,
        "genome_size_gb": 1.05,
        "gc_percent": 42.0,
        "morphology_type": "fowl",
        "sample_sequence": "ATGGGTCGTAAGCTGGATCTGCGCGAGTACGACCAGATCCTGGAGCGCATCGAGGAGATGCTTCAGCGCAAGGAGCTG",
        "description": "Galloanserae reference anchor, microchromosome architecture, high recombination."
    },
    {
        "id": "r_cucullatus",
        "name": "Raphus cucullatus",
        "common_name": "Dodo (Extinct ~1662)",
        "clade": "Aves (Columbiformes - Extinct)",
        "divergence_mya": 23.4,
        "genome_size_gb": 1.15,
        "gc_percent": 41.8,
        "morphology_type": "dodo",
        "sample_sequence": "ATGGGTCGCAAGCTGGACCTGCGCGAGTATGATCAGATCCTCGAGCGCATTGAGGAGATGCTGCAGCGTAAGGAGCTG",
        "description": "Flightless pigeon relative from Mauritius, de-extinction candidate."
    },
    {
        "id": "archaeopteryx",
        "name": "Archaeopteryx lithographica",
        "common_name": "Urvogel (Prehistoric ~150 Mya)",
        "clade": "Avialae (Jurassic Transition - Extinct)",
        "divergence_mya": 150.0,
        "genome_size_gb": 1.40,
        "gc_percent": 40.5,
        "morphology_type": "archaeopteryx",
        "sample_sequence": "ATGGGTCGGAAGCTGGACCTGCGGGAGTACGACCAGATCCTGGAGCGCATCGAGGAGATGCTGCAGCGCAAGGAGCTG",
        "description": "Transitional fossil between non-avian feathered dinosaurs and modern birds."
    }
]

class HybridizationRequest(BaseModel):
    species_a_id: str
    species_b_id: str
    hybrid_ratio: float = 0.5 # 0.0 = 100% A, 1.0 = 100% B
    divergence_mya: Optional[float] = 20.0
    custom_seq_a: Optional[str] = None
    custom_seq_b: Optional[str] = None

@app.get("/api/species-catalog")
def get_species_catalog():
    return {"catalog": SPECIES_CATALOG}

@app.post("/api/hybridize")
def hybridize_species(req: HybridizationRequest):
    """
    Simulates in-silico hybridization and ancestral genomic reconstruction (ASR).
    Synthesizes intermediate hybrid DNA and evaluates reading frame preservation & viability.
    """
    # Find species
    sp_a = next((s for s in SPECIES_CATALOG if s["id"] == req.species_a_id), SPECIES_CATALOG[0])
    sp_b = next((s for s in SPECIES_CATALOG if s["id"] == req.species_b_id), SPECIES_CATALOG[1])

    seq_a = req.custom_seq_a or sp_a["sample_sequence"]
    seq_b = req.custom_seq_b or sp_b["sample_sequence"]

    n = max(len(seq_a), len(seq_b))
    # Pad to equal length
    seq_a_pad = seq_a.ljust(n, 'N')
    seq_b_pad = seq_b.ljust(n, 'N')

    # Synthesize Hybrid Sequence using weighted probabilistic interpolation
    hybrid_chars = []
    ancestral_chars = []
    diff_count = 0

    ratio_b = max(0.0, min(1.0, float(req.hybrid_ratio)))
    ratio_a = 1.0 - ratio_b

    for i in range(n):
        b_a = seq_a_pad[i].upper()
        b_b = seq_b_pad[i].upper()

        if b_a != b_b and b_a != 'N' and b_b != 'N':
            diff_count += 1
            # Weighted hybrid selection
            chosen = b_b if (i % 100) / 100.0 < ratio_b else b_a
            # Ancestral consensus (MRCA) favors transition over transversion
            purines = {'A', 'G'}
            pyrimidines = {'C', 'T'}
            if (b_a in purines and b_b in purines) or (b_a in pyrimidines and b_b in pyrimidines):
                anc_chosen = b_a # Transition drift
            else:
                anc_chosen = b_a if ratio_a >= 0.5 else b_b
        else:
            chosen = b_a if b_a != 'N' else b_b
            anc_chosen = chosen

        hybrid_chars.append(chosen)
        ancestral_chars.append(anc_chosen)

    hybrid_seq = "".join(hybrid_chars)
    ancestral_seq = "".join(ancestral_chars)

    # Calculate GC content
    gc_count = sum(1 for c in hybrid_seq if c in ['G', 'C'])
    gc_pct = (gc_count / len(hybrid_seq) * 100.0) if hybrid_seq else 0.0

    # Translation to amino acids
    aa_list = []
    has_premature_stop = False
    for i in range(0, len(hybrid_seq) - 2, 3):
        codon = hybrid_seq[i:i+3]
        aa = GENETIC_CODE.get(codon, 'X')
        if aa == '*' and (i < len(hybrid_seq) - 6):
            has_premature_stop = True
        aa_list.append(aa)

    protein_seq = "".join(aa_list)

    # Compute biological viability score (0 - 100%)
    genetic_dist_pct = (diff_count / n * 100.0) if n > 0 else 0.0
    viability = 98.0 - (genetic_dist_pct * 0.4)
    if has_premature_stop:
        viability -= 35.0
    if not hybrid_seq.startswith("ATG"):
        viability -= 15.0
    viability_score = max(5.0, min(100.0, round(viability, 1)))

    return {
        "species_a": sp_a,
        "species_b": sp_b,
        "hybrid_ratio_a": round(ratio_a, 2),
        "hybrid_ratio_b": round(ratio_b, 2),
        "hybrid_sequence": hybrid_seq,
        "ancestral_mrca_sequence": ancestral_seq,
        "protein_translation": protein_seq,
        "length_bp": n,
        "genetic_distance_pct": round(genetic_dist_pct, 2),
        "gc_content_pct": round(gc_pct, 1),
        "viability_score": viability_score,
        "has_premature_stop": has_premature_stop,
        "morphology_mix": {
            "primary_morph": sp_a["morphology_type"] if ratio_a >= 0.5 else sp_b["morphology_type"],
            "secondary_morph": sp_b["morphology_type"] if ratio_a >= 0.5 else sp_a["morphology_type"],
            "blend_factor": ratio_b
        }
    }

# Mount static files
os.makedirs(STATIC_DIR, exist_ok=True)
app.mount("/static", StaticFiles(directory=STATIC_DIR), name="static")

@app.get("/")
def serve_index():
    return FileResponse(os.path.join(STATIC_DIR, "index.html"))

if __name__ == "__main__":
    import uvicorn
    print("\n[+] Starting PhyloGap Studio Web Dashboard on http://127.0.0.1:8085")
    uvicorn.run(app, host="127.0.0.1", port=8085)
