"""
PhyloGap - Gap Distribution & Ground Truth Audit
Parses DENTIST benchmark BED and JSON files to extract empirical gap statistics:
- Number of contigs and gaps
- Gap length distribution (min, 25%, median, 75%, max, mean)
- Gap length categories (<50bp, 50-200bp, 200-1000bp, >1000bp)
- Chromosome-level breakdown
"""

import os
import json
import numpy as np

BASE_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "benchmarks", "dentist_pilot")

def parse_bed_gaps(bed_path):
    """
    In DENTIST mapped-contigs.bed:
    Each line represents a mapped contig: chrom, start, end, name, score, strand ...
    Gaps are the inter-contig regions on the same chromosome/scaffold.
    """
    contigs_by_chrom = {}
    with open(bed_path, 'r') as f:
        for line in f:
            if not line.strip() or line.startswith('#') or line.startswith('track'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            if chrom not in contigs_by_chrom:
                contigs_by_chrom[chrom] = []
            contigs_by_chrom[chrom].append((start, end))

    # Calculate inter-contig gaps
    all_gaps = []
    gap_by_chrom = {}
    
    for chrom, intervals in contigs_by_chrom.items():
        intervals.sort(key=lambda x: x[0])
        gap_by_chrom[chrom] = []
        for i in range(len(intervals) - 1):
            prev_end = intervals[i][1]
            next_start = intervals[i+1][0]
            gap_len = next_start - prev_end
            if gap_len > 0:
                all_gaps.append(gap_len)
                gap_by_chrom[chrom].append(gap_len)

    return all_gaps, gap_by_chrom, sum(len(v) for v in contigs_by_chrom.values())

def analyze_gap_lengths(gaps):
    if not gaps:
        return {}
    arr = np.array(gaps)
    cat_less_50 = np.sum(arr < 50)
    cat_50_200 = np.sum((arr >= 50) & (arr < 200))
    cat_200_1000 = np.sum((arr >= 200) & (arr < 1000))
    cat_greater_1000 = np.sum(arr >= 1000)
    
    return {
        "count": len(arr),
        "min": int(np.min(arr)),
        "q25": float(np.percentile(arr, 25)),
        "median": float(np.median(arr)),
        "mean": float(np.mean(arr)),
        "q75": float(np.percentile(arr, 75)),
        "max": int(np.max(arr)),
        "categories": {
            "<50bp": (int(cat_less_50), float(cat_less_50 / len(arr) * 100)),
            "50-200bp": (int(cat_50_200), float(cat_50_200 / len(arr) * 100)),
            "200-1000bp": (int(cat_200_1000), float(cat_200_1000 / len(arr) * 100)),
            ">=1000bp": (int(cat_greater_1000), float(cat_greater_1000 / len(arr) * 100)),
        }
    }

def audit_json_report(json_path):
    with open(json_path, 'r') as f:
        data = json.load(f)
    return data

def main():
    print("=" * 60)
    print("PhyloGap - Gap Distribution & Ground Truth Audit")
    print("=" * 60)

    # 1. Drosophila melanogaster
    dmel_bed = os.path.join(BASE_DIR, "d_melanogaster", "assembly-reference.mapped-contigs.bed")
    dmel_json = os.path.join(BASE_DIR, "d_melanogaster", "closable-gaps.assembly-reference.reads-simulated-pb.json")
    
    if os.path.exists(dmel_bed):
        gaps, gap_by_chrom, n_contigs = parse_bed_gaps(dmel_bed)
        stats = analyze_gap_lengths(gaps)
        print(f"\n[1] Organism: Drosophila melanogaster (Fruit fly)")
        print(f"    - Mapped contigs: {n_contigs}")
        print(f"    - Inter-contig gaps: {stats['count']}")
        print(f"    - Gap length summary (bp): Min={stats['min']}, 25%={stats['q25']:.1f}, Median={stats['median']:.1f}, Mean={stats['mean']:.1f}, 75%={stats['q75']:.1f}, Max={stats['max']}")
        print(f"    - Stratification:")
        for cat, (cnt, pct) in stats['categories'].items():
            print(f"        * {cat:10s}: {cnt:5d} gaps ({pct:5.1f}%)")
        print(f"    - Top 5 scaffolds by gap count:")
        top_chroms = sorted(gap_by_chrom.items(), key=lambda x: len(x[1]), reverse=True)[:5]
        for c, g in top_chroms:
            print(f"        * {c}: {len(g)} gaps (median {np.median(g):.0f} bp)")

    if os.path.exists(dmel_json):
        jdata = audit_json_report(dmel_json)
        print(f"\n    [Closable Gaps Analysis - PacBio Simulated]:")
        print(f"    - JSON summary keys: {list(jdata.keys()) if isinstance(jdata, dict) else type(jdata)}")
        if isinstance(jdata, dict):
            for k in list(jdata.keys())[:6]:
                val = jdata[k]
                if isinstance(val, (int, float, str, bool)):
                    print(f"        * {k}: {val}")
                elif isinstance(val, list):
                    print(f"        * {k}: [list of {len(val)} items]")
                elif isinstance(val, dict):
                    print(f"        * {k}: [dict of {len(val)} keys]")

    # 2. Calypte anna (Avian)
    canna_bed = os.path.join(BASE_DIR, "c_anna", "assembly-reference.mapped-contigs.bed")
    if os.path.exists(canna_bed):
        gaps_c, gap_by_chrom_c, n_contigs_c = parse_bed_gaps(canna_bed)
        stats_c = analyze_gap_lengths(gaps_c)
        print(f"\n[2] Organism: Calypte anna (Anna's hummingbird - Avian focal)")
        print(f"    - Mapped contigs: {n_contigs_c}")
        print(f"    - Inter-contig gaps: {stats_c['count']}")
        print(f"    - Gap length summary (bp): Min={stats_c['min']}, 25%={stats_c['q25']:.1f}, Median={stats_c['median']:.1f}, Mean={stats_c['mean']:.1f}, 75%={stats_c['q75']:.1f}, Max={stats_c['max']}")
        print(f"    - Stratification:")
        for cat, (cnt, pct) in stats_c['categories'].items():
            print(f"        * {cat:10s}: {cnt:5d} gaps ({pct:5.1f}%)")
        print(f"    - Top 5 scaffolds by gap count:")
        top_chroms_c = sorted(gap_by_chrom_c.items(), key=lambda x: len(x[1]), reverse=True)[:5]
        for c, g in top_chroms_c:
            print(f"        * {c}: {len(g)} gaps (median {np.median(g):.0f} bp)")

    print("\n" + "=" * 60)

if __name__ == "__main__":
    main()
