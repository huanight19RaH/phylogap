"""
PhyloGap - Sequence-Level Exploratory Data Analysis (EDA)
Extracts exact ground truth gap sequences from assembly-reference.fasta.gz
using the benchmark coordinates in mapped-contigs.bed and closable-gaps.json:
1. Calculates global assembly statistics:
   - Scaffolds, total length, N50, GC content
2. Analyzes ground truth gap sequences:
   - GC-content of gap sequences vs GC-content of 100bp flanking regions
   - Homopolymer runs (poly-A/T/G/C >= 6bp)
   - Dinucleotide repeats (e.g. AT, GC, CA runs)
   - Gap length stratification vs biological complexity
"""

import os
import gzip
import re
import numpy as np

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "benchmarks", "dentist_pilot", "d_melanogaster")

def read_fasta_gz(filepath):
    header = None
    seq_parts = []
    with gzip.open(filepath, 'rt') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)
        if header is not None:
            yield header, "".join(seq_parts)

def compute_n50(lengths):
    lengths_sorted = sorted(lengths, reverse=True)
    total_len = sum(lengths_sorted)
    cumsum = 0
    for l in lengths_sorted:
        cumsum += l
        if cumsum >= total_len / 2.0:
            return l
    return 0

def calc_gc(seq):
    seq_upper = seq.upper()
    g = seq_upper.count('G')
    c = seq_upper.count('C')
    valid = len(seq) - seq_upper.count('N')
    if valid == 0:
        return 0.0
    return (g + c) / valid * 100.0

def detect_homopolymers(seq, min_run=6):
    pattern = re.compile(r'(A{' + str(min_run) + r',}|T{' + str(min_run) + r',}|C{' + str(min_run) + r',}|G{' + str(min_run) + r',})', re.IGNORECASE)
    matches = pattern.findall(seq)
    return len(matches) > 0, len(matches)

def detect_dinucleotide_repeats(seq, min_repeats=4):
    pattern = re.compile(r'((?:AT){' + str(min_repeats) + r',}|(?:TA){' + str(min_repeats) + r',}|(?:CG){' + str(min_repeats) + r',}|(?:GC){' + str(min_repeats) + r',}|(?:AC){' + str(min_repeats) + r',}|(?:TG){' + str(min_repeats) + r',})', re.IGNORECASE)
    matches = pattern.findall(seq)
    return len(matches) > 0, len(matches)

def main():
    ref_path = os.path.join(DATA_DIR, "assembly-reference.fasta.gz")
    bed_path = os.path.join(DATA_DIR, "assembly-reference.mapped-contigs.bed")

    if not os.path.exists(ref_path) or not os.path.exists(bed_path):
        print(f"[-] Required files not found in {DATA_DIR}.")
        return

    print("=" * 68)
    print("PhyloGap - Sequence-Level EDA & Biological Ground Truth Analysis")
    print("=" * 68)

    # 1. Load Reference Genomes
    print("\n[*] Loading assembly-reference.fasta.gz into memory...")
    ref_dict = {}
    ref_lengths = []
    total_ref_bp = 0
    for h, s in read_fasta_gz(ref_path):
        ref_dict[h] = s
        ref_lengths.append(len(s))
        total_ref_bp += len(s)
    
    print(f"    - Chromosomes/Scaffolds: {len(ref_dict)}")
    print(f"    - Total genome size: {total_ref_bp:,} bp ({total_ref_bp/(1024*1024):.2f} Mb)")
    print(f"    - Assembly N50: {compute_n50(ref_lengths):,} bp")
    all_chrom_seq = "".join(list(ref_dict.values()))
    genome_gc = calc_gc(all_chrom_seq)
    print(f"    - Baseline Genome GC Content: {genome_gc:.2f}%")

    # 2. Parse BED intervals and extract ground truth gaps
    print("\n[*] Parsing mapped-contigs.bed and extracting ground-truth gap sequences...")
    contigs_by_chrom = {}
    with open(bed_path, 'r') as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            if chrom not in contigs_by_chrom:
                contigs_by_chrom[chrom] = []
            contigs_by_chrom[chrom].append((start, end))

    gap_records = []
    for chrom, intervals in contigs_by_chrom.items():
        if chrom not in ref_dict:
            continue
        chrom_seq = ref_dict[chrom]
        intervals.sort(key=lambda x: x[0])
        for i in range(len(intervals) - 1):
            gap_start = intervals[i][1]
            gap_end = intervals[i+1][0]
            gap_len = gap_end - gap_start
            if gap_len > 0:
                true_seq = chrom_seq[gap_start:gap_end]
                left_flank = chrom_seq[max(0, gap_start - 100):gap_start]
                right_flank = chrom_seq[gap_end:min(len(chrom_seq), gap_end + 100)]
                flank_seq = left_flank + right_flank
                
                has_hp, hp_count = detect_homopolymers(true_seq, min_run=6)
                has_di, di_count = detect_dinucleotide_repeats(true_seq, min_repeats=4)
                
                gap_records.append({
                    "chrom": chrom,
                    "start": gap_start,
                    "end": gap_end,
                    "length": gap_len,
                    "true_seq": true_seq,
                    "gap_gc": calc_gc(true_seq),
                    "flank_gc": calc_gc(flank_seq),
                    "has_homopolymer": has_hp,
                    "hp_count": hp_count,
                    "has_dinucleotide": has_di,
                    "di_count": di_count,
                    "is_repetitive": has_hp or has_di
                })

    print(f"    - Extracted Ground-Truth Gaps: {len(gap_records):,}")

    # 3. Statistical Analysis of Ground Truth Sequences
    gap_gcs = [r["gap_gc"] for r in gap_records]
    flank_gcs = [r["flank_gc"] for r in gap_records]
    lengths = [r["length"] for r in gap_records]
    rep_flags = [r["is_repetitive"] for r in gap_records]
    hp_flags = [r["has_homopolymer"] for r in gap_records]
    di_flags = [r["has_dinucleotide"] for r in gap_records]

    print("\n" + "-" * 68)
    print("BIOLOGICAL FEATURES OF GROUND TRUTH GAP SEQUENCES:")
    print("-" * 68)
    print(f"1. GC Content Bias:")
    print(f"   - Whole Genome Baseline GC:       {genome_gc:.2f}%")
    print(f"   - Mean GC inside Ground-Truth Gaps: {np.mean(gap_gcs):.2f}% (Std: {np.std(gap_gcs):.2f}%)")
    print(f"   - Mean GC in 100bp Flanking Regions: {np.mean(flank_gcs):.2f}% (Std: {np.std(flank_gcs):.2f}%)")
    delta_gc = np.array(gap_gcs) - np.array(flank_gcs)
    print(f"   - Mean Delta GC (Gap - Flank):      {np.mean(delta_gc):+.2f}%")
    
    extreme_gc = np.sum((np.array(gap_gcs) < 25) | (np.array(gap_gcs) > 60))
    print(f"   - Gaps with Extreme GC (<25% or >60%): {extreme_gc} / {len(gap_records)} ({extreme_gc/len(gap_records)*100:.1f}%)")

    print(f"\n2. Low-Complexity & Repetitive Sequence Enrichment:")
    print(f"   - Gaps with Homopolymer runs (>=6bp):   {sum(hp_flags)} / {len(gap_records)} ({sum(hp_flags)/len(gap_records)*100:.1f}%)")
    print(f"   - Gaps with Dinucleotide repeats (>=4x): {sum(di_flags)} / {len(gap_records)} ({sum(di_flags)/len(gap_records)*100:.1f}%)")
    print(f"   - Total Repetitive / Low-Complexity:     {sum(rep_flags)} / {len(gap_records)} ({sum(rep_flags)/len(gap_records)*100:.1f}%)")

    print(f"\n3. Gap Length Stratification Breakdown:")
    for min_l, max_l, label in [(0, 50, "<50bp"), (50, 200, "50-200bp"), (200, 1000, "200-1000bp"), (1000, 100000, ">=1000bp")]:
        sub = [r for r in gap_records if min_l <= r["length"] < max_l]
        if sub:
            sub_gc = np.mean([r["gap_gc"] for r in sub])
            sub_rep = sum(r["is_repetitive"] for r in sub) / len(sub) * 100
            print(f"   - {label:10s} (N={len(sub):4d}): Mean GC={sub_gc:5.1f}%, Repetitive Enrichment={sub_rep:5.1f}%")

    print("\n" + "=" * 68)

if __name__ == "__main__":
    main()
