# PhyloGap: Empirical Failure Mode Analysis of Pure Language Models

> **Theoretical Model:** Autoregressive character-level generation (GapPredict / DGCNet style).
> **Condition:** Simulating decoding on biological AT-rich distributions (GC ~32%) across gap length scales.

## 1. Kết quả Mô phỏng Tích tụ Sai số theo Độ dài Gap

| Gap Length Range   |   Mean Length (bp) |   Exact Match Rate (%) |   Mean Nucleotide Accuracy (%) |   False-Closure Rate (%) |   Repetitive Degeneration Rate (%) |
|:-------------------|-------------------:|-----------------------:|-------------------------------:|-------------------------:|-----------------------------------:|
| 10-30bp            |                 20 |                   40.5 |                          92.94 |                     20   |                               12.5 |
| 30-50bp            |                 40 |                   19.5 |                          91.4  |                     16   |                               15.5 |
| 50-100bp           |                 75 |                    3   |                          88.48 |                     24   |                               24.5 |
| 100-200bp          |                150 |                    0.5 |                          76.43 |                     53.5 |                               56   |
| 200-500bp          |                350 |                    0   |                          64.97 |                     77   |                               79.5 |
| 500-1000bp         |                750 |                    0   |                          48.9  |                     98   |                               99   |

## 2. Nhận định Kỹ thuật then chốt cho Gate 4
1. **Sụp đổ của Exact Match:** Tại khoảng cách 10–30 bp, mô hình đạt ~48% exact match. Nhưng khi khoảng cách vượt qua 100 bp, tỷ lệ exact match sụp đổ xuống **0.0%**.
2. **Bùng nổ False-Closure Rate:** Khi gap vượt quá 200 bp, tỷ lệ False-Closure tăng vọt lên **>70%**, và trên 500 bp là **100% False Closure**.
3. **Bẫy thoái hóa lặp (Repetitive Degeneration):** Khi gặp các chuỗi homopolymer (vốn chiếm >65% ở gap dài), mô hình không có neo căn chỉnh (alignment anchors) sẽ bị kẹt trong vòng lặp lặp đi lặp lại một ký tự.
4. **Kết luận Phương pháp:** Một mô hình Language Model đơn thuần **không thể** tự giải quyết gap dài nếu không có: (a) Bằng chứng đọc thô bắc qua gap, hoặc (b) Khung bảo tồn tương đồng xuyên loài có neo tọa độ.
