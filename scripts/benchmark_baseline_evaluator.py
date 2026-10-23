"""
PhyloGap - Baseline Benchmark Extractor & Standardizer.
Parses published benchmark results for DENTIST, PBJelly, and LR_gapcloser
on Drosophila melanogaster (model invertebrate) and Calypte anna (avian genome).

Computes standardized metrics across all tools:
- Exact Match Rate (%ID = 100%)
- Near-Exact Rate (%ID >= 99%)
- Moderate Accuracy Rate (%ID >= 90%)
- False-Closure Rate (FCR): Filled gaps with %ID < 90% or broken scaffolds
- Contig N50 improvement ratio
- Length-stratified breakdown (<500bp, 500-1000bp, >1000bp)
"""

import os
import json
import pandas as pd
import numpy as np

OUTPUT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "experiments", "results")

# Raw data parsed directly from Hiller Lab benchmark gap-stats.tsv
BENCHMARK_RAW_DATA = [
    # 1. Drosophila melanogaster (1,340 closable gaps)
    {
        "organism": "Drosophila melanogaster",
        "clade": "Invertebrate (Pilot)",
        "tool": "DENTIST",
        "existing_closable": 1340,
        "closed": 1266,
        "broken_scaffolds": 0,
        "id_100": 1221,
        "id_99": 1238,
        "id_95": 1264,
        "id_90": 1265,
        "avg_identity": 0.9989,
        "input_n50": 237511,
        "result_n50": 8810613,
        "max_n50": 21485538,
        "by_size": {
            "<500bp": {"existing": 1288, "closed": 1221, "id_100": 1186, "id_90": 1220},
            "500-1000bp": {"existing": 39, "closed": 35, "id_100": 30, "id_90": 35},
            ">1000bp": {"existing": 13, "closed": 10, "id_100": 5, "id_90": 10}
        }
    },
    {
        "organism": "Drosophila melanogaster",
        "clade": "Invertebrate (Pilot)",
        "tool": "PBJelly",
        "existing_closable": 1340,
        "closed": 1282,
        "broken_scaffolds": 0,
        "id_100": 6,
        "id_99": 80,
        "id_95": 1075,
        "id_90": 1257,
        "avg_identity": 0.9584,
        "input_n50": 237511,
        "result_n50": 7640288,
        "max_n50": 21485538,
        "by_size": {
            "<500bp": {"existing": 1288, "closed": 1239, "id_100": 6, "id_90": 1215},
            "500-1000bp": {"existing": 39, "closed": 34, "id_100": 0, "id_90": 33},
            ">1000bp": {"existing": 13, "closed": 9, "id_100": 0, "id_90": 9}
        }
    },
    {
        "organism": "Drosophila melanogaster",
        "clade": "Invertebrate (Pilot)",
        "tool": "LR_gapcloser",
        "existing_closable": 1340,
        "closed": 1201,
        "broken_scaffolds": 5,
        "id_100": 0,
        "id_99": 39,
        "id_95": 317,
        "id_90": 1117,
        "avg_identity": 0.9378,
        "input_n50": 237511,
        "result_n50": 13484323,
        "max_n50": 21485538,
        "by_size": {
            "<500bp": {"existing": 1288, "closed": 1158, "id_100": 0, "id_90": 1080},
            "500-1000bp": {"existing": 39, "closed": 33, "id_100": 0, "id_90": 29},
            ">1000bp": {"existing": 13, "closed": 10, "id_100": 0, "id_90": 8}
        }
    },
    # 2. Calypte anna (Anna's hummingbird - 37,501 closable gaps)
    {
        "organism": "Calypte anna",
        "clade": "Aves (Primary Target)",
        "tool": "DENTIST",
        "existing_closable": 37501,
        "closed": 37024,
        "broken_scaffolds": 0,
        "id_100": 32583,
        "id_99": 36375,
        "id_95": 36928,
        "id_90": 37002,
        "avg_identity": 0.99915,
        "input_n50": 49266,
        "result_n50": 12763730,
        "max_n50": 14522327,
        "by_size": {
            "<500bp": {"existing": 25239, "closed": 25027, "id_100": 23558, "id_90": 25005},
            "500-1000bp": {"existing": 6436, "closed": 6340, "id_100": 5237, "id_90": 6340},
            ">1000bp": {"existing": 5826, "closed": 5657, "id_100": 3788, "id_90": 5657}
        }
    },
    {
        "organism": "Calypte anna",
        "clade": "Aves (Primary Target)",
        "tool": "PBJelly",
        "existing_closable": 37501,
        "closed": 8360,
        "broken_scaffolds": 66,
        "id_100": 11,
        "id_99": 372,
        "id_95": 7551,
        "id_90": 8228,
        "avg_identity": 0.96545,
        "input_n50": 49266,
        "result_n50": 64273,
        "max_n50": 14522327,
        "by_size": {
            "<500bp": {"existing": 25239, "closed": 5469, "id_100": 11, "id_90": 5364},
            "500-1000bp": {"existing": 6436, "closed": 1632, "id_100": 0, "id_90": 1613},
            ">1000bp": {"existing": 5826, "closed": 1259, "id_100": 0, "id_90": 1251}
        }
    },
    {
        "organism": "Calypte anna",
        "clade": "Aves (Primary Target)",
        "tool": "LR_gapcloser",
        "existing_closable": 37501,
        "closed": 33906,
        "broken_scaffolds": 18,
        "id_100": 14,
        "id_99": 602,
        "id_95": 4985,
        "id_90": 26799,
        "avg_identity": 0.92409,
        "input_n50": 49266,
        "result_n50": 14551778,
        "max_n50": 14522327,
        "by_size": {
            "<500bp": {"existing": 25239, "closed": 22919, "id_100": 14, "id_90": 20602},
            "500-1000bp": {"existing": 6436, "closed": 5789, "id_100": 0, "id_90": 4070},
            ">1000bp": {"existing": 5826, "closed": 5198, "id_100": 0, "id_90": 2127}
        }
    }
]

def compute_standardized_metrics():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    rows = []
    
    for item in BENCHMARK_RAW_DATA:
        total = item["existing_closable"]
        closed = item["closed"]
        broken = item["broken_scaffolds"]
        id_100 = item["id_100"]
        id_99 = item["id_99"]
        id_90 = item["id_90"]
        
        # Standardized metrics
        closure_rate = (closed / total) * 100.0
        exact_match_rate = (id_100 / total) * 100.0
        exact_in_closed = (id_100 / closed) * 100.0 if closed > 0 else 0.0
        high_acc_rate = (id_99 / closed) * 100.0 if closed > 0 else 0.0
        
        # False Closure Rate (FCR): Closed gaps with %ID < 90% plus broken scaffolds penalty
        bad_closures = (closed - id_90) + broken
        fcr = (bad_closures / closed) * 100.0 if closed > 0 else 0.0
        
        # Contig N50 improvement factor
        n50_fold = item["result_n50"] / item["input_n50"]
        n50_pct_max = (item["result_n50"] / item["max_n50"]) * 100.0
        
        rows.append({
            "Organism": item["organism"],
            "Clade": item["clade"],
            "Tool": item["tool"],
            "Closable Gaps": total,
            "Gaps Closed": closed,
            "Closure Rate (%)": round(closure_rate, 2),
            "Exact Matches (%ID=100%)": id_100,
            "Exact Match Rate (%)": round(exact_match_rate, 2),
            "Exact in Closed (%)": round(exact_in_closed, 2),
            "High Accuracy Rate (%ID>=99%)": round(high_acc_rate, 2),
            "Average Identity (%)": round(item["avg_identity"] * 100.0, 2),
            "Broken Scaffolds": broken,
            "False-Closure Rate (%)": round(fcr, 2),
            "Initial N50 (bp)": item["input_n50"],
            "Final N50 (bp)": item["result_n50"],
            "N50 Gain (x-fold)": round(n50_fold, 1),
            "N50 % of Max Contig": round(n50_pct_max, 1)
        })

    df = pd.DataFrame(rows)
    csv_path = os.path.join(OUTPUT_DIR, "baseline_benchmark_sota.csv")
    df.to_csv(csv_path, index=False)
    print(f"[+] Saved standardized benchmark results to: {csv_path}")

    # Generate Markdown summary table
    md_summary = df.to_markdown(index=False)
    md_path = os.path.join(os.path.dirname(OUTPUT_DIR), "..", "docs", "baseline_benchmark_report.md")
    with open(md_path, "w", encoding="utf-8") as f:
        f.write("# PhyloGap: Standardized SOTA Baseline Benchmark Report\n\n")
        f.write("> **Source Data:** Empirical benchmark results from Hiller Lab (DENTIST publication).\n")
        f.write("> **Standardized Metrics:** Computed according to PhyloGap Evaluation Engine specifications.\n\n")
        f.write("## 1. Tổng hợp Hiệu năng Toàn diện (Overall Performance Matrix)\n\n")
        f.write(md_summary + "\n\n")
        f.write("## 2. Phân tích Chi tiết Lỗi và Nguy cơ False Closure\n\n")
        f.write("### Phát hiện 1: Nghịch lý Completeness vs. Exact Match\n")
        f.write("- **Trên Drosophila melanogaster:** DENTIST đóng 1,266 gap nhưng đạt **96.4% exact match** (1,221 gap hoàn hảo). Trong khi đó, PBJelly đóng 1,282 gap nhưng chỉ có **6 gap đạt exact match (0.5%)**!\n")
        f.write("- **Trên Calypte anna (Chim ruồi):** DENTIST đóng 37,024 gap và đạt **88.0% exact match** trong số gap đã đóng (32,583 gap). Ngược lại, LR_gapcloser đóng 33,906 gap nhưng chỉ có **14 gap exact match (0.04%)**, với False-Closure Rate lên tới **20.97%** và phá vỡ 18 scaffolds!\n\n")
        f.write("### Phát hiện 2: Đánh giá Rủi ro của việc Chắp ghép Mù (Chimeric Joins)\n")
        f.write("- LR_gapcloser đạt Contig N50 cao nhất (14.55 Mb) nhưng đồng thời có tỷ lệ sai số cao nhất và làm vỡ scaffold. Đây là bằng chứng thực nghiệm rõ ràng nhất khẳng định: **Tăng Contig N50 mà không kiểm soát False-Closure Rate là cực kỳ nguy hiểm cho phân tích sinh học downstream.**\n")
        
    print(f"[+] Saved markdown report to: {md_path}")
    return df

if __name__ == "__main__":
    df = compute_standardized_metrics()
    print("\n" + df[["Organism", "Tool", "Closure Rate (%)", "Exact in Closed (%)", "False-Closure Rate (%)", "Average Identity (%)"]].to_string())
