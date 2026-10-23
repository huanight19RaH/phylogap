# PhyloGap: Standardized SOTA Baseline Benchmark Report

> **Source Data:** Empirical benchmark results from Hiller Lab (DENTIST publication).
> **Standardized Metrics:** Computed according to PhyloGap Evaluation Engine specifications.

## 1. Tổng hợp Hiệu năng Toàn diện (Overall Performance Matrix)

| Organism                | Clade                 | Tool         |   Closable Gaps |   Gaps Closed |   Closure Rate (%) |   Exact Matches (%ID=100%) |   Exact Match Rate (%) |   Exact in Closed (%) |   High Accuracy Rate (%ID>=99%) |   Average Identity (%) |   Broken Scaffolds |   False-Closure Rate (%) |   Initial N50 (bp) |   Final N50 (bp) |   N50 Gain (x-fold) |   N50 % of Max Contig |
|:------------------------|:----------------------|:-------------|----------------:|--------------:|-------------------:|---------------------------:|-----------------------:|----------------------:|--------------------------------:|-----------------------:|-------------------:|-------------------------:|-------------------:|-----------------:|--------------------:|----------------------:|
| Drosophila melanogaster | Invertebrate (Pilot)  | DENTIST      |            1340 |          1266 |              94.48 |                       1221 |                  91.12 |                 96.45 |                           97.79 |                  99.89 |                  0 |                     0.08 |             237511 |          8810613 |                37.1 |                  41   |
| Drosophila melanogaster | Invertebrate (Pilot)  | PBJelly      |            1340 |          1282 |              95.67 |                          6 |                   0.45 |                  0.47 |                            6.24 |                  95.84 |                  0 |                     1.95 |             237511 |          7640288 |                32.2 |                  35.6 |
| Drosophila melanogaster | Invertebrate (Pilot)  | LR_gapcloser |            1340 |          1201 |              89.63 |                          0 |                   0    |                  0    |                            3.25 |                  93.78 |                  5 |                     7.41 |             237511 |         13484323 |                56.8 |                  62.8 |
| Calypte anna            | Aves (Primary Target) | DENTIST      |           37501 |         37024 |              98.73 |                      32583 |                  86.89 |                 88.01 |                           98.25 |                  99.91 |                  0 |                     0.06 |              49266 |         12763730 |               259.1 |                  87.9 |
| Calypte anna            | Aves (Primary Target) | PBJelly      |           37501 |          8360 |              22.29 |                         11 |                   0.03 |                  0.13 |                            4.45 |                  96.55 |                 66 |                     2.37 |              49266 |            64273 |                 1.3 |                   0.4 |
| Calypte anna            | Aves (Primary Target) | LR_gapcloser |           37501 |         33906 |              90.41 |                         14 |                   0.04 |                  0.04 |                            1.78 |                  92.41 |                 18 |                    21.01 |              49266 |         14551778 |               295.4 |                 100.2 |

## 2. Phân tích Chi tiết Lỗi và Nguy cơ False Closure

### Phát hiện 1: Nghịch lý Completeness vs. Exact Match
- **Trên Drosophila melanogaster:** DENTIST đóng 1,266 gap nhưng đạt **96.4% exact match** (1,221 gap hoàn hảo). Trong khi đó, PBJelly đóng 1,282 gap nhưng chỉ có **6 gap đạt exact match (0.5%)**!
- **Trên Calypte anna (Chim ruồi):** DENTIST đóng 37,024 gap và đạt **88.0% exact match** trong số gap đã đóng (32,583 gap). Ngược lại, LR_gapcloser đóng 33,906 gap nhưng chỉ có **14 gap exact match (0.04%)**, với False-Closure Rate lên tới **20.97%** và phá vỡ 18 scaffolds!

### Phát hiện 2: Đánh giá Rủi ro của việc Chắp ghép Mù (Chimeric Joins)
- LR_gapcloser đạt Contig N50 cao nhất (14.55 Mb) nhưng đồng thời có tỷ lệ sai số cao nhất và làm vỡ scaffold. Đây là bằng chứng thực nghiệm rõ ràng nhất khẳng định: **Tăng Contig N50 mà không kiểm soát False-Closure Rate là cực kỳ nguy hiểm cho phân tích sinh học downstream.**
