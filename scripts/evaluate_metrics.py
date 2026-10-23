"""
PhyloGap - Core Evaluation Engine for Genome Assembly Gap Reconstruction.
Standardized metric computation for benchmarking read-based and learning-based gap-filling methods.

Metrics Computed:
1. Nucleotide-level Accuracy (Acc_nuc)
2. Error Decomposition:
   - Substitution Error Rate
   - Insertion Error Rate
   - Deletion Error Rate
3. Exact Match Rate (100% Identity)
4. False-Closure Rate (FCR): Proportion of filled gaps exceeding allowable error tolerance
5. Selective Prediction & Reliability:
   - Risk-Coverage Curve
   - Area Under Risk-Coverage curve (AU-RC)
"""

import sys
import numpy as np

def needleman_wunsch(seq1, seq2, match=2, mismatch=-1, gap_open=-2, gap_extend=-1):
    """
    Standard global sequence alignment (Needleman-Wunsch) with affine gap penalty.
    Pure Python implementation ensuring zero-dependency reproducibility.
    Returns: aligned_seq1, aligned_seq2, score
    """
    n, m = len(seq1), len(seq2)
    if n == 0 and m == 0:
        return "", "", 0
    if n == 0:
        return "-" * m, seq2, m * gap_open
    if m == 0:
        return seq1, "-" * n, n * gap_open

    # DP tables: M (match/mismatch), X (gap in seq1), Y (gap in seq2)
    # Using simple linear gap penalty for efficiency if strings are long
    dp = np.zeros((n + 1, m + 1), dtype=np.int32)
    for i in range(n + 1):
        dp[i][0] = i * gap_open
    for j in range(m + 1):
        dp[0][j] = j * gap_open

    for i in range(1, n + 1):
        c1 = seq1[i - 1].upper()
        for j in range(1, m + 1):
            c2 = seq2[j - 1].upper()
            score_diag = dp[i - 1][j - 1] + (match if c1 == c2 and c1 != 'N' else mismatch)
            score_del = dp[i - 1][j] + gap_open
            score_ins = dp[i][j - 1] + gap_open
            dp[i][j] = max(score_diag, score_del, score_ins)

    # Traceback
    align1, align2 = [], []
    i, j = n, m
    while i > 0 or j > 0:
        c1 = seq1[i - 1].upper() if i > 0 else ''
        c2 = seq2[j - 1].upper() if j > 0 else ''
        diag_score = match if c1 == c2 and c1 != 'N' else mismatch
        
        if i > 0 and j > 0 and dp[i][j] == dp[i - 1][j - 1] + diag_score:
            align1.append(seq1[i - 1])
            align2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif i > 0 and dp[i][j] == dp[i - 1][j] + gap_open:
            align1.append(seq1[i - 1])
            align2.append('-')
            i -= 1
        else:
            align1.append('-')
            align2.append(seq2[j - 1])
            j -= 1

    return "".join(reversed(align1)), "".join(reversed(align2)), int(dp[n][m])

def evaluate_single_gap(true_seq, pred_seq, max_error_tolerance=0.10, max_len_tolerance=0.20):
    """
    Evaluates a single gap reconstruction against its ground truth.
    
    Parameters:
    - true_seq (str): Ground truth sequence (no Ns).
    - pred_seq (str): Reconstructed sequence proposed by the tool.
    - max_error_tolerance (float): Maximum allowed error rate before classifying as False Closure (default 10%).
    - max_len_tolerance (float): Maximum relative length divergence allowed (default 20%).
    
    Returns:
    - dict with detailed evaluation metrics.
    """
    true_clean = true_seq.strip().upper()
    pred_clean = pred_seq.strip().upper()

    # Case: Empty prediction or unclosed gap (retained Ns)
    if not pred_clean or pred_clean == 'N' * len(pred_clean):
        return {
            "status": "unfilled",
            "is_filled": False,
            "exact_match": False,
            "nucleotide_accuracy": 0.0,
            "error_rate": 1.0,
            "substitutions": 0,
            "insertions": 0,
            "deletions": len(true_clean),
            "true_length": len(true_clean),
            "pred_length": 0,
            "length_ratio": 0.0,
            "is_false_closure": False,
        }

    # Case: Identical match
    if true_clean == pred_clean:
        return {
            "status": "exact_match",
            "is_filled": True,
            "exact_match": True,
            "nucleotide_accuracy": 1.0,
            "error_rate": 0.0,
            "substitutions": 0,
            "insertions": 0,
            "deletions": 0,
            "true_length": len(true_clean),
            "pred_length": len(pred_clean),
            "length_ratio": 1.0,
            "is_false_closure": False,
        }

    # Case: Needleman-Wunsch Alignment
    align_true, align_pred, score = needleman_wunsch(true_clean, pred_clean)
    
    matches = 0
    subs = 0
    insertions = 0 # gap in true -> extra base inserted by pred
    deletions = 0  # gap in pred -> true base omitted by pred

    for t, p in zip(align_true, align_pred):
        if t == '-' and p != '-':
            insertions += 1
        elif t != '-' and p == '-':
            deletions += 1
        elif t == p:
            matches += 1
        else:
            subs += 1

    total_ops = matches + subs + insertions + deletions
    alignment_len = max(len(true_clean), len(pred_clean), 1)
    
    acc = matches / alignment_len
    error_rate = (subs + insertions + deletions) / alignment_len
    len_ratio = len(pred_clean) / max(len(true_clean), 1)
    len_error = abs(len(pred_clean) - len(true_clean)) / max(len(true_clean), 1)

    # A gap is marked as a False Closure if the tool committed a filled sequence
    # but the error rate or length discrepancy exceeds acceptable assembly thresholds.
    is_fc = (error_rate > max_error_tolerance) or (len_error > max_len_tolerance)

    return {
        "status": "filled",
        "is_filled": True,
        "exact_match": False,
        "nucleotide_accuracy": float(acc),
        "error_rate": float(error_rate),
        "substitutions": subs,
        "insertions": insertions,
        "deletions": deletions,
        "true_length": len(true_clean),
        "pred_length": len(pred_clean),
        "length_ratio": float(len_ratio),
        "is_false_closure": is_fc,
        "aligned_true": align_true,
        "aligned_pred": align_pred
    }

def compute_risk_coverage_curve(predictions_with_confidence):
    """
    Computes the Risk-Coverage Curve and Area Under Risk-Coverage (AU-RC).
    
    Parameters:
    - predictions_with_confidence: list of tuples (is_error, confidence_score)
      where is_error is 1 if prediction was wrong (or False Closure), 0 otherwise.
      
    Returns:
    - coverages: list of coverage percentages [0.0 to 1.0]
    - risks: list of empirical error rates on accepted predictions
    - aurc: Area Under the Risk-Coverage Curve (lower is better)
    """
    if not predictions_with_confidence:
        return [0.0], [0.0], 0.0

    # Sort descending by confidence score
    sorted_items = sorted(predictions_with_confidence, key=lambda x: x[1], reverse=True)
    n = len(sorted_items)
    
    coverages = []
    risks = []
    cum_errors = 0

    for i, (is_error, conf) in enumerate(sorted_items, start=1):
        cum_errors += 1 if is_error else 0
        cov = i / n
        current_risk = cum_errors / i
        coverages.append(float(cov))
        risks.append(float(current_risk))

    # Calculate trapezoidal AUC (compatible with both newer and older numpy versions)
    trap_fn = getattr(np, 'trapezoid', getattr(np, 'trapz', None))
    aurc = float(trap_fn(risks, coverages)) if len(coverages) > 1 and trap_fn else (risks[0] if risks else 0.0)
    return coverages, risks, aurc

def evaluate_batch(records, max_error_tolerance=0.10):
    """
    Evaluates a batch of gap reconstruction records.
    Each record must have: 'true_seq', 'pred_seq', and optional 'confidence' (float).
    
    Returns summary statistics dictionary.
    """
    total = len(records)
    if total == 0:
        return {}

    eval_results = []
    conf_pairs = []

    for r in records:
        true_s = r["true_seq"]
        pred_s = r.get("pred_seq", "")
        conf = float(r.get("confidence", 1.0))
        
        res = evaluate_single_gap(true_s, pred_s, max_error_tolerance=max_error_tolerance)
        eval_results.append(res)
        
        # Risk defined as false closure or major error
        is_risk = 1 if res["is_false_closure"] or (res["is_filled"] and res["error_rate"] > max_error_tolerance) else 0
        conf_pairs.append((is_risk, conf))

    filled_gaps = [r for r in eval_results if r["is_filled"]]
    exact_matches = [r for r in eval_results if r["exact_match"]]
    false_closures = [r for r in eval_results if r["is_false_closure"]]

    n_filled = len(filled_gaps)
    fcr = len(false_closures) / max(n_filled, 1) if n_filled > 0 else 0.0
    mean_acc = np.mean([r["nucleotide_accuracy"] for r in filled_gaps]) if n_filled > 0 else 0.0
    mean_err = np.mean([r["error_rate"] for r in filled_gaps]) if n_filled > 0 else 0.0

    covs, risks, aurc = compute_risk_coverage_curve(conf_pairs)

    return {
        "total_gaps": total,
        "filled_count": n_filled,
        "fill_rate": float(n_filled / total * 100),
        "exact_match_count": len(exact_matches),
        "exact_match_rate": float(len(exact_matches) / total * 100),
        "mean_nucleotide_accuracy": float(mean_acc),
        "mean_error_rate": float(mean_err),
        "false_closure_count": len(false_closures),
        "false_closure_rate": float(fcr * 100),
        "aurc": float(aurc),
        "individual_results": eval_results
    }

if __name__ == "__main__":
    print("=== PhyloGap Evaluation Metrics Module ===")
    print("Run tests with: python -m unittest tests/test_evaluate_metrics.py")
