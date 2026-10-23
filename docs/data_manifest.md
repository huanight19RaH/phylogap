# PhyloGap: Dataset Manifest & Ground Truth Specification (v0.1)

> **Document Status:** Active Audit · **Date:** 23/10/2026  
> **Repository Path:** `docs/data_manifest.md`  
> **Data Root:** `data/benchmarks/`

---

## 1. Tổng quan Kiến trúc Dữ liệu Đa tầng (Multi-Tier Data Design)

Dự án PhyloGap xây dựng ma trận dữ liệu 3 tầng (Tier 1: Pilot & Tái lập; Tier 2: Benchmark Chuẩn hóa Đa loài; Tier 3: Miền Chim B10K):

```mermaid
graph TD
    T1[Tier 1: Pilot Baselines] --> T2[Tier 2: Multi-Species Ground Truth]
    T2 --> T3[Tier 3: Avian Generalization & Evolutionary Decay]
    
    T1_1[D. melanogaster<br/>1,374 gaps | dm6] --> T1
    T1_2[A. thaliana<br/>Col-0 reference] --> T1
    
    T2_1[C. anna Hummingbird<br/>38,996 realistic gaps] --> T2
    T2_2[H. sapiens GRCh38<br/>Controlled / Real reads] --> T2
    
    T3_1[B10K PRJNA545868<br/>Sister species clade] --> T3
    T3_2[B10K 363-bird WGA<br/>Phylogenetic distances] --> T3
```

---

## 2. Kết quả Thẩm định Thực nghiệm Bộ Dữ liệu Pilot (Empirical Audit)

Dựa trên dữ liệu tọa độ contig và ground truth gap tải về từ máy chủ MPI-CBG Hiller Lab (`assembly-reference.mapped-contigs.bed` và `closable-gaps.*.json`):

### 2.1. Bảng Thống kê So sánh Cơ cấu Gap: *Drosophila melanogaster* vs. *Calypte anna*

| Chỉ số Thống kê | *Drosophila melanogaster* (Pilot) | *Calypte anna* (Anna's Hummingbird - Avian) | Ý nghĩa Kỹ thuật đối với PhyloGap |
| :--- | :---: | :---: | :--- |
| **Số lượng Mapped Contigs** | 3,244 | 39,155 | Hệ gen chim bị phân mảnh cao gấp >10 lần |
| **Số lượng Gap thực tế ($N$-gaps)** | **1,374** | **38,996** | Quy mô benchmark chim lớn gấp gần 30 lần |
| **Độ dài Gap nhỏ nhất (Min)** | 11 bp | 10 bp | Cả hai đều có gap vi mô do lỗi căn chỉnh |
| **Phân vị 25% ($Q_1$)** | 33.0 bp | 49.0 bp | 75% gap ở chim dài hơn 49 bp |
| **Trung vị Độ dài Gap (Median)** | **45.0 bp** | **208.0 bp** | **Gap ở chim dài gấp >4.6 lần so với ruồi giấm** |
| **Trung bình Độ dài Gap (Mean)** | 84.7 bp | 498.2 bp | Gap chim có phân phối đuôi dài (*heavy-tailed*) |
| **Phân vị 75% ($Q_3$)** | 86.0 bp | 619.0 bp | 25% gap ở chim dài vượt quá 619 bp |
| **Độ dài Gap lớn nhất (Max)** | 1,104 bp | 59,392 bp | Tồn tại các siêu gap do lặp hoặc dị nhiễm sắc |

### 2.2. Phân tầng Chi tiết theo Kích thước Gap (Length Stratification)

```text
Drosophila melanogaster (1,374 gaps)
├── <50 bp      :  733 gaps (53.3%)  ███████████████████████████
├── 50–200 bp   :  529 gaps (38.5%)  ███████████████████
├── 200–1000 bp :  109 gaps ( 7.9%)  ████
└── >=1000 bp   :    3 gaps ( 0.2%)  ▏

Calypte anna (Anna's hummingbird) (38,996 gaps)
├── <50 bp      : 9,845 gaps (25.2%)  █████████████
├── 50–200 bp   : 9,382 gaps (24.1%)  ████████████
├── 200–1000 bp : 14,606 gaps (37.5%) ███████████████████
└── >=1000 bp   : 5,163 gaps (13.2%)  ███████
```

> **Phát hiện quan trọng cho Thiết kế Thuật toán:**  
> Ở sinh vật mô hình nhỏ (*Drosophila*), **91.8%** số gap có kích thước dưới 200 bp. Các mô hình Deep Learning trước đây như **GapPredict (2021)** hay **DLGapCloser (2024)** chỉ thử nghiệm trên nấm và sinh vật nhỏ nên dễ đạt độ tương đồng cao.  
> Tuy nhiên, ở chim (*Calypte anna*), hơn **50.7%** số gap dài hơn 200 bp (trong đó 13.2% gap dài $\ge 1000$ bp). Các thuật toán dùng Beam Search cục bộ sẽ gặp bùng nổ không gian trạng thái (*exponential search space*) và sụt giảm độ chính xác nghiêm trọng.

---

## 3. Cấu trúc Ground Truth & Thông tin Định danh (Ground Truth Spec)

Trong tệp JSON `closable-gaps.assembly-reference.reads-simulated-pb.json`:
Mỗi gap được định nghĩa chính xác theo cấu trúc đối tượng:

```json
{
  "fromContig": 1,
  "toContig": 2,
  "gapSize": 36,
  "mappedInterval": {
    "tag": 1,
    "begin": 187143,
    "end": 187179,
    "contigId": 1
  },
  "spanningReads": [4532, 13436, 15411, 19565, 26711, 29557, 43988, 64899]
}
```

* `fromContig`, `toContig`: ID của 2 contig liền kề trong test assembly.
* `gapSize`: Khoảng cách nucleotide chưa xác định (độ dài chuỗi `N`).
* `mappedInterval`: Tọa độ chính xác `[begin, end)` trên reference genome hoàn chỉnh (Ground Truth sequence).
* `spanningReads`: Danh sách các reads thực sự bắc qua gap trong tập dữ liệu giải trình tự (dùng để kiểm chứng Track 1: Read-supported).

---

## 4. Kế hoạch Thu thập Dữ liệu B10K (Avian Clade Selection & Curated Accessions)

Để trả lời **RQ3 (Ảnh hưởng của khoảng cách tiến hóa)**, PhyloGap chọn 4 cấp độ phân tách tiến hóa từ dự án B10K (PRJNA545868) và VGP, được thẩm định và xác thực tự động qua NCBI Datasets API (chi tiết xem tại [avian_clade_accessions.md](file:///d:/THStudy/UniversityStudy/Research/phylogap/docs/avian_clade_accessions.md) và tệp manifest có cấu trúc [avian_b10k_manifest.csv](file:///d:/THStudy/UniversityStudy/Research/phylogap/data/manifests/avian_b10k_manifest.csv)):

| Phân tầng (Tier) | Loài Mục tiêu (Target) & Accession | Loài Tham chiếu (Reference) & Accession | Phân tách (Mya) | Cấp độ Lắp ráp | Mục đích Thực nghiệm RQ3 |
| :--- | :--- | :--- | :---: | :---: | :--- |
| **Tier 1 (Sister Species)** | *Lonchura striata*<br/>`GCF_046129695.1` (bLonStr1.mat) | *Lonchura oryzivora*<br/>`GCA_056776775.1` (ASM5677677v1) | **~5.2** | Chromosome | Độ tương đồng cực cao; kiểm chứng trần hiệu năng của homology-guided filling |
| **Tier 2 (Same Family)** | *Taeniopygia guttata*<br/>`GCF_003957565.2` (bTaeGut1.4.pri) | *Lonchura striata*<br/>`GCF_046129695.1` (bLonStr1.mat) | **~11.5** | Chromosome | Đánh giá sai số do substitution & vi đột biến indel trong cùng họ Estrildidae |
| **Tier 2b (Alt Family Ref)** | *Taeniopygia guttata*<br/>`GCF_003957565.2` (bTaeGut1.4.pri) | *Poephila acuticauda*<br/>`GCA_057534655.1` (Pacu_hecki_V1) | **~14.1** | Chromosome | Nhóm đối chứng bổ trợ độc lập cho nhánh chim sẻ cỏ Úc |
| **Tier 3 (Same Order)** | *Gallus gallus*<br/>`GCF_016699485.2` (GRCg7b) | *Numida meleagris*<br/>`GCF_002078875.1` (NumMel1.0) | **~46.2** | Chromosome | Khảo sát mức độ suy thoái căn chỉnh tương đồng giữa các họ trong bộ Galliformes |
| **Tier 4 (Distant Orders)** | *Calypte anna*<br/>`GCF_003957555.1` (bCalAnn1_v1.p) | *Gallus gallus*<br/>`GCF_016699485.2` (GRCg7b) | **~88.0** | Chromosome | Xác định điểm gãy (*breakdown point*) khi loài tham chiếu quá xa (Galloanserae vs Neoaves) |

> **Ghi chú Thẩm định Kỹ thuật:** *Taeniopygia bichenovii* ban đầu được thay thế bởi *Lonchura oryzivora* do *T. bichenovii* chưa có assembly de novo trên NCBI. Tất cả các assembly được chọn đều là dữ liệu chính thức từ VGP/B10K, đạt cấp độ Chromosome-scale với Contig N50 từ 8.9 Mb đến 36.1 Mb.

