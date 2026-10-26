# PhyloGap: Kế Hoạch Thực Nghiệm Chi Tiết Pha 5 (Phase 5 Detailed Experimental Plan)

> **Dự án:** PhyloGap — Tái cấu trúc khoảng trống bộ gene nhận biết bằng chứng (Evidence-Aware Genome Gap Reconstruction)  
> **Tài liệu:** Kế hoạch thực nghiệm mở rộng & Thử nghiệm đối đầu SOTA (Phase 5 Experimental Plan)  
> **Tác giả:** Nhóm Nghiên cứu PhyloGap · Người thực hiện: Thái Hưng Nguyễn  
> **Ngày lập:** 2026-10-26 · **Phạm vi:** In Silico Computational Genomics (Dry lab strictly)  
> **Nguyên tắc:** Sử dụng 100% dữ liệu thật, kiểm định chặt chẽ, zero emoji, khả năng tái lập 100%.

---

## 1. Mục Tiêu Thực Nghiệm (Scientific Objectives)

Chuyển giao từ mô phỏng lý thuyết (RQ1–RQ5 đã nghiệm thu) sang thực nghiệm định lượng trên dữ liệu hệ gene thực tế, giải quyết 3 bài toán lớn:
1. **Bài toán Đánh đổi Completeness vs Correctness:** Chứng minh rằng các công cụ SOTA đọc sâu (DENTIST) bị tê liệt (tỷ lệ bỏ trống gap >85%) khi reads bị thưa (1–2x), trong khi mô hình ngôn ngữ (GapPredict) và homology thô gây ra lỗi đóng giả (False Closures >90%).
2. **Khảo sát Hiệu năng Pareto của PhyloGap:** Chứng minh thuật toán Evidence Fusion kết hợp Selective Gating vừa khôi phục được gap bằng reads thưa vừa giữ False-Closure Rate < 3.0%.
3. **Hiệu Chuẩn Tốc Độ Đột Biến Thích Ứng (Adaptive Rate Lambda):** Chứng minh tính vượt trội của việc cá nhân hóa $\lambda$ theo vùng chức năng sinh học (coding exon, intron, intergenic) và độ lệch GC so với scalar $\lambda$ toàn cục.

---

## 2. Nguồn Dữ Liệu Thực Nghiệm (Authentic Real-World Datasets)

Tuyệt đối không sử dụng dữ liệu tự tạo ngẫu nhiên thiếu căn cứ sinh học. Dữ liệu thực nghiệm bao gồm:

### 2.1. Dữ liệu Chuẩn Drosophila melanogaster (DENTIST Benchmark)
* **Tệp dữ liệu:**
  * `data/benchmarks/dentist_pilot/d_melanogaster/assembly-reference.fasta.gz` (~40.2 MB): Bộ gene chuẩn chromosome-level (ground truth).
  * `data/benchmarks/dentist_pilot/d_melanogaster/closable-gaps.assembly-reference.reads-simulated-pb.json`: Tọa độ gap thực tế, flanking sequences và chỉ số reads bắc qua.
* **Mục đích:** Đánh giá trên các vị trí gap sinh học thực tế có phân phối GC và độ phức tạp cấu trúc (repeats, homopolymers) phản ánh đúng tự nhiên.

### 2.2. Dữ liệu Hệ Gene Chim B10K (Bird 10,000 Genomes Project / VGP)
Trích xuất theo manifest chuẩn `data/manifests/avian_b10k_manifest.csv`:
* **Tầng 1 (Sister species, ~5.2 Mya):**
  * Target: *Lonchura striata* (GCF_046129695.1, PacBio Revio HiFi, Chromosome-level).
  * Reference: *Lonchura oryzivora* (GCA_056776775.1, PacBio HiFi, Chromosome-level).
* **Tầng 2 (Same family, ~11.5 Mya):**
  * Target: *Taeniopygia guttata* (GCF_003957565.2, VGP Phase 1, Chromosome-level).
  * Reference: *Lonchura striata* (GCF_046129695.1).
* **Tầng 3 (Same order, ~46.2 Mya):**
  * Target: *Gallus gallus* (GCF_016699485.2, GRCg7b, Chromosome-level).
  * Reference: *Numida meleagris* (GCF_002078875.1, INRA, Chromosome-level).
* **Tầng 4 (Distant orders, ~88.0 Mya):**
  * Target: *Calypte anna* (GCF_003957555.1, VGP Phase 1, Chromosome-level).
  * Reference: *Gallus gallus* (GCF_016699485.2).

---

## 3. Thiết Kế Ma Trận Thực Nghiệm (Experimental Matrix)

Thực hiện tổng cộng **720 thử nghiệm có kiểm soát** (Controlled Gap Trials):

```
                                  [720 THỬ NGHIỆM ĐỐI ĐẦU]
                                             │
             ┌───────────────────────────────┼───────────────────────────────┐
             │                               │                               │
       [4 TẦNG TIẾN HÓA]             [3 KỊCH BẢN READS]              [3 ĐỘ DÀI GAP]
    • Tầng 1: 5.2 Mya (Chi)        • Đầy đủ (>20x Depth)           • Ngắn: 50 bp
    • Tầng 2: 11.5 Mya (Họ)        • Thưa: 1-2x (Non-model)        • Vừa: 100 bp
    • Tầng 3: 46.2 Mya (Bộ)        • Không có reads (0x)           • Dài: 200 bp
    • Tầng 4: 88.0 Mya (Lớp)
                                             │
             └───────────────────────────────┼───────────────────────────────┘
                                             │
                              [4 TRƯỜNG PHÁI PHƯƠNG PHÁP LUẬN]
                  1. DENTIST (Read-Based Consensus & Repeat Filter)
                  2. GapPredict (DNA Language Model Autoregressive)
                  3. Raw Homology (Direct Orthologous Scaffold Transfer)
                  4. PhyloGap Engine (Evidence Fusion + Selective Gating)
```

---

## 4. Các Chỉ Số Đo Lường & Tiêu Chí Nghiệm Thu (Evaluation Framework)

1. **Closure Rate (%):** Tỷ lệ phần trăm gap được công cụ chấp nhận điền.
2. **Exact-Match Rate (%):** Tỷ lệ phần trăm gap được khôi phục chính xác 100% so với Ground Truth.
3. **Nucleotide Accuracy ($Acc_{nuc}$ %):** Tỷ lệ nucleotide khớp chính xác sau căn chỉnh toàn cục Needleman-Wunsch.
4. **False-Closure Rate ($FCR$ %):**
   $$FCR = \frac{\text{Số gap đã điền nhưng bị sai lệch vượt ngưỡng cho phép}}{\text{Tổng số gap được công cụ chấp nhận điền}} \times 100\%$$
   *Tiêu chí an toàn: FCR < 5.0% là an toàn, FCR > 50.0% là mức nguy hiểm (chimeric assembly).*
5. **Selective Gating Rejection Rate (%):** Tỷ lệ gap được bảo vệ (giữ nguyên ký tự `N`) khi khoảng cách tiến hóa vượt mốc $t_{crit}$.

---

## 5. Kế Hoạch Triển Khai Từng Bước (Execution Roadmap)

* **Bước 1: Kiểm định nền tảng (Base Invariants):**
  * Chạy test suite `tests/test_adaptive_and_head_to_head.py` và `tests/test_evidence_fusion.py` để xác thực toàn vẹn toán học.
* **Bước 2: Chạy Thực Nghiệm Đối Đầu (Head-to-Head Execution):**
  * Chạy script `scripts/benchmark_head_to_head_sota.py` trên 720 trials.
  * Xuất dữ liệu thô ra `experiments/results/head_to_head_sota_trials.csv`.
  * Xuất bảng tổng hợp ra `experiments/results/head_to_head_sota_summary.csv`.
  * Tạo báo cáo khoa học tại `docs/head_to_head_sota_report.md`.
* **Bước 3: Tích hợp Cấu hình Điều phối Đa Tài khoản Kaggle:**
  * Cập nhật `kaggle_runner_ui/data/projects.yaml` với 3 job Phase 5 cho 3 tài khoản (`thung192`, `TThanh13`, `RaH1111`).
* **Bước 4: Cập nhật Sổ tay Nghiên cứu & Đồng bộ Notion:**
  * Ghi nhận kết quả thực nghiệm vào `RESEARCH_PLAN.md`.
  * Đồng bộ commit hash và chỉ số thực nghiệm lên Notion Project Portfolio và Weekly Action Board.
