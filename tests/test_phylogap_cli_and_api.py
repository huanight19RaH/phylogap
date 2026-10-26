"""
Unit tests for PhyloGap CLI pipeline and FastAPI Web Server endpoints.
"""

import os
import sys
import pytest

# Ensure repo root is on sys.path
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from scripts.phylogap_cli import reconstruct_gap_pipeline, DEMO_PRESETS
from phylogap_ui.server import app
from fastapi.testclient import TestClient

client = TestClient(app)


def test_demo_presets_structure():
    """Verify built-in demo presets are properly structured."""
    assert 1 in DEMO_PRESETS
    assert 2 in DEMO_PRESETS
    assert 3 in DEMO_PRESETS

    for tier_id, preset in DEMO_PRESETS.items():
        assert "tier" in preset
        assert "target_species" in preset
        assert "reference_species" in preset
        assert "flank_left" in preset
        assert "flank_right" in preset
        assert "gap_truth" in preset
        assert "divergence_mya" in preset
        assert "sparse_reads" in preset
        assert len(preset["flank_left"]) >= 20
        assert len(preset["flank_right"]) >= 20


def test_reconstruct_gap_pipeline_tier2():
    """Test full CLI pipeline reconstruction on Zebra Finch preset."""
    p2 = DEMO_PRESETS[2]
    # Build synthetic reference locus with evolved sequence
    from scripts.cross_species_homology_engine import AvianEvolutionarySimulator
    sim = AvianEvolutionarySimulator(seed=42)
    target = p2["flank_left"] + p2["gap_truth"] + p2["flank_right"]
    evolved, _ = sim.evolve_sequence(target, divergence_time_mya=p2["divergence_mya"])
    ref_locus = ("NNNN" * 5) + evolved + ("NNNN" * 5)

    res = reconstruct_gap_pipeline(
        flank_left=p2["flank_left"],
        flank_right=p2["flank_right"],
        homology_reference=ref_locus,
        divergence_mya=p2["divergence_mya"],
        expected_gap_len=len(p2["gap_truth"]),
        sparse_reads=p2["sparse_reads"],
        ground_truth_gap=p2["gap_truth"]
    )

    assert "is_accepted" in res
    assert "final_sequence" in res
    assert "confidence_score" in res
    assert res["is_accepted"] is True
    assert res["ground_truth_eval"] is not None
    # Accuracy should be high (>90%) with read correction
    assert res["ground_truth_eval"]["nucleotide_accuracy"] > 0.90
    assert res["read_corrections_applied"] > 0


def test_api_health_endpoint():
    """Verify /api/health endpoint returns 200 and valid metadata."""
    resp = client.get("/api/health")
    assert resp.status_code == 200
    data = resp.json()
    assert data["status"] == "healthy"
    assert data["project"] == "PhyloGap"
    assert data["t_crit_mya"] == 13.2


def test_api_presets_endpoint():
    """Verify /api/presets returns 3 validated presets."""
    resp = client.get("/api/presets")
    assert resp.status_code == 200
    data = resp.json()
    assert "presets" in data
    assert len(data["presets"]) == 3
    tiers = [p["tier_id"] for p in data["presets"]]
    assert tiers == [1, 2, 3]


def test_api_risk_calc_endpoint():
    """Verify /api/risk-calc accurately reflects evolutionary decay and danger zone."""
    # Test safe tier
    resp_safe = client.get("/api/risk-calc?mya=4.0")
    assert resp_safe.status_code == 200
    data_safe = resp_safe.json()
    assert data_safe["risk_tier"] == "SAFE_CONGENERIC"
    assert data_safe["is_exceeding_t_crit"] is False
    assert data_safe["expected_nucleotide_identity"] > 95.0

    # Test danger tier
    resp_danger = client.get("/api/risk-calc?mya=46.2")
    assert resp_danger.status_code == 200
    data_danger = resp_danger.json()
    assert data_danger["risk_tier"] == "DANGER_ORDER_CLASS"
    assert data_danger["is_exceeding_t_crit"] is True
    assert data_danger["homology_weight_w_hom"] == 0.0


def test_api_reconstruct_endpoint():
    """Verify /api/reconstruct returns complete provenance ribbon and metrics."""
    p1 = DEMO_PRESETS[1]
    from scripts.cross_species_homology_engine import AvianEvolutionarySimulator
    sim = AvianEvolutionarySimulator(seed=101)
    target = p1["flank_left"] + p1["gap_truth"] + p1["flank_right"]
    evolved, _ = sim.evolve_sequence(target, divergence_time_mya=p1["divergence_mya"])
    ref_locus = ("NNNN" * 3) + evolved + ("NNNN" * 3)

    payload = {
        "flank_left": p1["flank_left"],
        "flank_right": p1["flank_right"],
        "homology_reference": ref_locus,
        "divergence_mya": p1["divergence_mya"],
        "expected_gap_len": len(p1["gap_truth"]),
        "sparse_reads": p1["sparse_reads"],
        "ground_truth_gap": p1["gap_truth"]
    }

    resp = client.post("/api/reconstruct", json=payload)
    assert resp.status_code == 200
    data = resp.json()
    assert "provenance" in data
    assert len(data["provenance"]) > 0
    assert "is_accepted" in data
    assert data["is_accepted"] is True
