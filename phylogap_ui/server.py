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
