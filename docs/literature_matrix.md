# PhyloGap: Comprehensive Literature Matrix & SOTA Benchmark Synthesis

> **Version:** 1.0 · **Date:** 23/10/2026  
> **Status:** Active Reference for Literature Review & Novelty Validation  
> **Repository Path:** `docs/literature_matrix.md`

---

## 1. Executive Summary & Research Positioning (Định vị Nghiên cứu)

Bài toán **Genome Assembly Gap Filling** đã trải qua 3 làn sóng công nghệ:
1. **Làn sóng 1 (2010–2018): De Bruijn Graph & Short Reads** (*Sealer, GapCloser, GapFiller*). Ưu điểm: tiết kiệm bộ nhớ; Nhược điểm: hoàn toàn bất lực trước các vùng lặp (*repetitive elements*) và khoảng trống dài hơn chiều dài read ($>150-300$ bp).
2. **Làn sóng 2 (2019–2023): Long Reads & Consensus Verification** (*DENTIST, PBJelly, LR_GapCloser, TGS-GapCloser*). DENTIST đã tạo nên bước ngoặt khi ưu tiên **độ chính xác tuyệt đối** thông qua căn chỉnh nhận biết vùng lặp và cơ chế xác minh chống ghép sai (*false closure*). Tuy nhiên, phương pháp phụ thuộc chặt chẽ vào độ phủ giải trình tự cao ($>20\times$) của long reads chất lượng tốt.
3. **Làn sóng 3 (2021–2026): Deep Learning, Homology & Heuristic Search** (*GapPredict, DLGapCloser, DL-GapFilling*). Nhóm nghiên cứu của Yu Chen & Gang Wang (2024, 2026) đã đưa CNN-BiLSTM, tìm kiếm *Wave-Beam* và gần đây là *BeamStar* cùng cơ chế *PredictionFilter* vào giải quyết gap trên genome nấm và thực vật.

### Khoảng trống Nghiên cứu mà PhyloGap Khai thác (The PhyloGap Research Gap)
Mặc dù các mô hình AI năm 2024–2026 phát triển rất nhanh, khảo sát kỹ thuật chỉ ra **4 hạn chế cốt tử chưa được giải quyết trong y văn**:
* **Hạn chế 1 (Thiếu thẩm định tiến hóa xuyên loài):** Cả *DLGapCloser (2024)* và *DL-GapFilling (2026)* đều chỉ kiểm thử trên vi nấm hoặc thực vật đơn lẻ; chưa từng có công trình nào định lượng mối tương quan giữa sự suy giảm độ chính xác tái dựng với **khoảng cách phát sinh loài (phylogenetic distance)** của genome tham chiếu.
* **Hạn chế 2 (Khoảng cách Gap ở Động vật có xương sống):** Thực nghiệm của chúng tôi trên *Calypte anna* (chim ruồi) cho thấy **50.7%** gap dài trên 200 bp (và 13.2% gap $\ge 1000$ bp). Các mô hình ngôn ngữ nucleotide dự đoán từng base (*autoregressive / beam search*) sẽ bị tích tụ sai số (*compounding error*) theo cấp số nhân khi khoảng cách vượt quá 200 bp.
* **Hạn chế 3 (Ảo giác vs. False Closure):** Chưa có nghiên cứu nào phân tích đường cong đánh đổi rủi ro (*Risk–Coverage Curve*) để cho phép mô hình AI **từ chối đưa ra dự đoán (selective prediction / abstention)** khi bằng chứng đọc và bằng chứng tương đồng mâu thuẫn nhau.
* **Hạn chế 4 (Thiếu Ground Truth chuẩn hóa):** Nhiều bài báo AI tự đánh giá bằng tỷ lệ % gap được điền trên draft assembly mà không đối chiếu với ground truth có kiểm soát tại các tọa độ gap thực tế như DENTIST đã làm.

---

## 2. Bảng Ma trận So sánh Toàn diện (Master SOTA Comparison Matrix)

| Tiêu chí So sánh | **DENTIST (2022)** | **GapPredict (2021)** | **DLGapCloser (2024)** | **DL-GapFilling (2026)** | **GapSuite (2026)** | **Comprehensive Eval (2024)** | **PhyloGap (Đề xuất)** |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Tác giả & Đơn vị** | Ludwig et al. (MPI-CBG, Hiller Lab) | Li et al. (BCGSC, Birol Lab) | Chen, Wang, Zhang (NEFU) | Chen, Wang, Wang, Wang (NEFU) | Liang et al. | Zhao, Liu, Pan | Nghiên cứu này |
| **Tạp chí & Định danh** | *GigaScience*, giab100 | *IEEE/ACM TCBB*, PMC8772386 | *IJMS*, 25(15):8502 (PMC11313336) | *Briefings in Bioinformatics*, bbag007 | *Bioinformatics / Tools* | *Genes*, PMC10815858 | Computational Plan v0.1 |
| **Mô hình / Thuật toán Lõi** | Long-read alignment, Repeat masking, DAG consensus, Validation filter | Character-level LSTM Language Model | DGCNet (1D-CNN + BiLSTM) + Sealer re-assembly | DFillingNet (Deep Filling NN) + Structural Priors | Gap-Aid (extension) + Gap-Graph (graph path) | Benchmark so sánh 7 công cụ long-read | Evidence-aware Fusion + Selective Prediction |
| **Thuật toán Tìm kiếm Chuỗi** | Pairwise local alignment (Dazzler) | Greedy / Argmax nucleotide probability | **Wave-Beam Search** (Expand-Shrink levels) | **BeamStar Search (BSCEA)** + Cost Function | Graph traversal (De Bruijn / Assembly graph) | Theo từng công cụ (MUMmer, BLAT, Graph) | Calibrated Beam Search / Rejection Threshold |
| **Dữ liệu Đầu vào Bắt buộc** | Assembly draft + Long reads (PacBio/ONT) | Assembly draft + Flanking regions + Reads | Flanking regions + Reads + Homologous Genome | Flanking regions + Structural sequence features | Raw reads + Assembly graph / Contigs | Tùy từng công cụ khảo sát | Assembly draft + Multi-scale Homology + Reads |
| **Ground Truth Thiết kế** | **Realistic Gap Coordinates** map sang Reference chuẩn | In silico masking trên draft | Homology alignment consistency qua Exonerate | Benchmark assembly chuẩn | T2T complete reference genomes | Reference assemblies (Haploid/Tetraploid) | **Realistic Gap Mapping (DENTIST-style) + B10K** |
| **Sinh vật & Quy mô Đã Đánh giá** | *Drosophila*, *Arabidopsis*, Chim ruồi (*C. anna*), Người | Vi khuẩn, Nấm men (*S. cerevisiae*) | 4 genome nhỏ: *S. cerevisiae*, *S. pombe*, *N. crassa*, *M. pusilla* | Hệ gen thực vật và vi tảo (*Arabidopsis*, Algae) | Thực vật (*Arabidopsis*, Lúa), Người | Ruồi giấm, Nấm men, Đa bội (*Arabidopsis*) | **Hệ gen Chim (B10K: Passeriformes, Galliformes, Apodiformes)** |
| **Chỉ số Đánh giá Cốt lõi** | Exact closure, Correct insertions, False-closure rate, QUAST | Sequence similarity, Gap closure % | Gap closed %, N count reduction, QUAST | Number of filled gaps, Flanking consistency | T2T contiguity, Completeness, Manual curation reduction | Completeness vs. Accuracy, Ploidy impact | **Nucleotide Acc, False-Closure Rate (FCR), Risk-Coverage Curve** |
| **Failure Modes Đã Nhận diện** | Cần read depth cao ($>15\times$); bất lực khi không có reads | Dễ sinh chuỗi ngẫu nhiên; không giải quyết được gap $>100$ bp | Tỷ lệ mismatch/indel tăng cao ở vùng lặp; chỉ chạy được genome nhỏ | Phức tạp trong huấn luyện; chưa kiểm chứng tính tổng quát xuyên loài | Phụ thuộc vào chất lượng assembly graph ban đầu | Sai khác lớn giữa đơn bội và đa bội; FGAP vượt trội | Suy thoái dự đoán khi khoảng cách tiến hóa tăng |

---

## 3. Phân tích Kỹ thuật Chi tiết Từng Công trình Hạt nhân

### 3.1. DENTIST (Ludwig et al., 2022) — Chuẩn mực về Độ Chính xác & Thiết kế Benchmark
* **Điểm đột phá:**
  1. Thay vì che ngẫu nhiên (*random masking*), DENTIST lấy vị trí gap thực tế từ draft assembly bị phân mảnh, sau đó đối chiếu tọa độ sang reference genome hoàn chỉnh. Điều này giữ nguyên các đặc tính sinh học tự nhiên của gap (vùng lặp, GC content cực đoan).
  2. Áp dụng cơ chế **Validation Filter**: Mỗi đoạn chèn (*insert*) được kiểm tra chéo với các reads hỗ trợ; nếu phát hiện bằng chứng mâu thuẫn (chimeric joins hoặc reads bắc qua không đồng thuận), DENTIST thà giữ nguyên `N` chứ không điền sai.
* **Hạn chế đối với PhyloGap:** DENTIST không sử dụng thông tin bảo tồn tương đồng xuyên loài (*homology*) và không có mô hình ngôn ngữ sinh học để nội suy khi reads bị thiếu hụt hoặc độ sâu giải trình tự cục bộ bằng 0.

### 3.2. DLGapCloser (Chen et al., 2024) — Deep Learning kết hợp Homology
* **Kiến trúc DGCNet:**
  * Lớp 1D-CNN (32 channels, kernel size 3, ReLU) để trích xuất đặc trưng cục bộ k-mer.
  * Lớp Bi-LSTM (512 chiều ẩn) để nắm bắt ngữ cảnh 2 chiều từ hai phía sườn gap (*flanking sequences*).
  * Lớp Dense phân loại 4 nucleotide (A, C, G, T).
* **Thuật toán Wave-Beam Search:**
  * Giải quyết nhược điểm cắt tỉa sớm (*early pruning*) của Beam Search truyền thống bằng cách chia chu kỳ:
    * *Expansion phase* (Levels 1–3, 5–6): Giữ lại toàn bộ các nhánh ứng viên.
    * *Contraction phase* (Levels 4, 7): Cắt tỉa dựa trên xác suất tích lũy (*cumulative probability sum*).
* **Hạn chế đối với PhyloGap:**
  1. Chỉ thử nghiệm trên 4 hệ gen kích thước rất nhỏ ($<40$ Mb).
  2. Bước cuối cùng lại ném toàn bộ dự đoán ngược trở lại vào công cụ *Sealer* (vốn là de Bruijn graph của short reads), làm mất đi tính độc lập của mô hình và giới hạn khả năng điền các gap dài.
  3. Mismatch và indel rate được chính tác giả ghi nhận là tăng cao rõ rệt ở các vùng lặp.

### 3.3. DL-GapFilling (Chen et al., 2026) — Cải tiến Tối ưu Tìm kiếm & Lọc Ứng viên
* **Điểm mới:**
  * Nâng cấp lên mạng *DFillingNet* và thuật toán *BeamStar (BSCEA)* kết hợp các tiên nghiệm cấu trúc bộ gen (*genomic structural priors*).
  * Bổ sung mô-đun *PredictionFilter* để loại bỏ các đoạn dự đoán có độ tự tin thấp trước khi đưa vào assembly.
* **Hạn chế đối với PhyloGap:**
  1. Nghiên cứu hoàn toàn khu trú trong hệ gen thực vật (*plant genomes*), nơi cấu trúc đa bội (*polyploidy*) và chuyển vị gen có đặc thù rất khác so với hệ gen chim (vốn có kích thước nhỏ gọn $\sim 1-1.2$ Gb nhưng mật độ lặp cục bộ rất phức tạp).
  2. Chưa khảo sát kịch bản **Cross-Species Evolutionary Divergence Decay** (khi khoảng cách tiến hóa thay đổi từ 5 Mya đến 80 Mya).

### 3.4. GapSuite (2026) — Hướng tiếp cận T2T Thực dụng
* **Cơ chế:** Kết hợp *Gap-Aid* (mở rộng đầu đọc thô) và *Gap-Graph* (tìm đường đi trên đồ thị assembly) nhằm hướng tới Telomere-to-Telomere (T2T).
* **Giá trị tham khảo:** Cung cấp giải pháp xử lý đồ thị lắp ráp thay vì chỉ coi chuỗi DNA như văn bản NLP thuần túy.

---

## 4. Bảng Đối chiếu Điểm Mới (Novelty Justification for PhyloGap)

Để bảo đảm vượt qua **Gate 1 (Tuần 2)** và **Gate 4 (Tuần 8)**, bảng dưới đây làm rõ sự khác biệt bản chất giữa PhyloGap và các công trình SOTA:

| Trục Nghiên cứu | SOTA Hiện hành (DL-GapFilling, DLGapCloser, DENTIST) | Hướng Đột phá của PhyloGap |
| :--- | :--- | :--- |
| **Miền Sinh học** | Nấm, Vi tảo, Thực vật, Sinh vật mô hình (*Drosophila*, Người). | **Miền Chim (*Aves* - B10K):** Hệ gen có tỷ lệ synteny cao nhưng chứa các vi nhiễm sắc thể (*microchromosomes*) giàu GC và gap dài. |
| **Biến số Khoảng cách Tiến hóa** | Bỏ qua hoặc chỉ dùng 1 genome tương đồng gần nhất không xác định niên đại. | **Phân tầng Tiến hóa 4 Cấp độ:** Khảo sát có hệ thống sự suy thoái độ chính xác từ Sister species ($\sim 3$ Mya) đến Distant orders ($>80$ Mya). |
| **Phân phối Kích thước Gap** | Tối ưu cho gap ngắn ($<100-200$ bp). | **Thiết kế chuyên biệt cho gap dài:** Giải quyết dải gap 200–1000 bp và $>1000$ bp (chiếm $>50.7\%$ gap ở chim). |
| **Xử lý Mâu thuẫn Bằng chứng** | Chọn chuỗi có xác suất ngôn ngữ cao nhất hoặc ghép mù vào Sealer. | **Evidence-Aware Fusion:** Trọng số hóa bằng chứng giữa Raw Reads, Homology Alignment và Language Model Priors. |
| **Độ an toàn của Assembly** | Tối đa hóa số lượng gap được điền (Completeness-first). | **Safety-first (Minimizing False-Closure Rate):** Thiết lập cơ chế Selective Prediction với đường cong kiểm soát sai số Risk–Coverage. |

---

## 5. Kết luận cho Giai đoạn Tiếp theo

1. **Về mặt Lý thuyết:** Y văn đã chứng minh Deep Learning có thể tham gia vào gap filling, nhưng các giải pháp hiện tại đều "mù tiến hóa" (*evolutionary-agnostic*) và thất bại trên các gap dài ở động vật bậc cao.
2. **Về mặt Thực nghiệm:** Việc tải và kiểm tra dữ liệu *Calypte anna* đã xác nhận giả thuyết rằng gap ở chim có đặc tính phân phối hoàn toàn khác biệt so với sinh vật mô hình nhỏ.
3. **Mục tiêu ưu tiên Tuần 3–4:** Tiến hành xây dựng script đánh giá độc lập (`scripts/evaluate_metrics.py`) để đo đạc nucleotide accuracy, false closures và sẵn sàng chạy thử nghiệm baseline đầu tiên trên *D. melanogaster*.
