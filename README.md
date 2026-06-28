# Nghiên cứu kỹ thuật ước lượng kênh DACE hỗ trợ bởi công suất đỉnh kết hợp giảm PAPR trong hệ thống MIMO-OFDM

> **Đồ án Thiết kế I** — Trường Điện - Điện tử, Đại học Bách khoa Hà Nội  
> **Sinh viên:** Trần Văn Thắng — Lớp Điện tử 01 - K68  
> **Giảng viên hướng dẫn:** CN. Ma Việt Đức  
> **Năm học:** 2025–2026

---

## Giới thiệu

Đồ án này nghiên cứu và xây dựng sơ đồ **ước lượng kênh hỗ trợ dữ liệu (DACE)** kết hợp kỹ thuật giảm tỉ số công suất đỉnh trên công suất trung bình (**PAPR**) trong hệ thống **MIMO-OFDM**, sử dụng kỹ thuật nén **Modified Gamma Correction Companding (MGCC)**.

**Ý tưởng cốt lõi:** tái sử dụng các subcarrier có công suất đỉnh cao — nguyên nhân gây PAPR cao — làm pilot bổ sung tại phía thu để cải thiện ước lượng kênh, không cần tăng mật độ pilot hay phản hồi từ phía thu.

---

## Cấu trúc thư mục

```
├── Test_1.m                  # Hàm mô phỏng chính (Monte Carlo)
├── mgcc_compress.m            # Hàm nén MGCC
├── mgcc_expand.m              # Hàm giãn MGCC (nghịch đảo đại số đúng)
├── Plot_PAPR.m                # Vẽ CCDF PAPR (Fig 4)
├── Plot_test_1.m              # Vẽ MSE SISO LS BPSK L=16 (Fig 5)
├── Plot_test_1_full.m         # Vẽ toàn bộ Fig 6–15
└── README.md
```

---

## Yêu cầu

- **MATLAB** R2020b trở lên
- **Communications Toolbox**
- **Signal Processing Toolbox**

---

## Cách chạy

Đặt tất cả file `.m` vào cùng một thư mục, sau đó chạy theo thứ tự:

```matlab
% Bước 1: Vẽ CCDF PAPR (Fig 4)
run('Plot_PAPR.m')

% Bước 2: Vẽ MSE SISO LS BPSK L=16 (Fig 5)
run('Plot_test_1.m')

% Bước 3: Vẽ toàn bộ Fig 6-15
run('Plot_test_1_full.m')
```

Hoặc gọi trực tiếp hàm mô phỏng để lấy số liệu một điểm SNR:

```matlab
% SISO, Rayleigh, LS, L=16, BPSK, rho=0
[mse_base, mse_tx, mse_rx, ber] = Test_1(10, 'Rayleigh', 'LS', 16, 2, 0);

% SIMO 1x2, Rician, LS, L=16, 4QAM, rho=0
[mse_base, mse_tx, mse_rx, ber] = Test_1(10, 'Rician', 'LS', 16, 4, 0, 1, 2);
```

---

## Tham số hệ thống

| Tham số | Giá trị |
|---------|---------|
| Số subcarrier N | 256 |
| Số pilot Np | 16 (comb-type) |
| Số subcarrier đỉnh Nrsc | 16 |
| Số tap kênh L | 4, 8, 16 |
| Mô hình kênh | Rayleigh, Rician (K=10) |
| Điều chế | BPSK, 4QAM, 8PSK |
| Tham số MGCC — Rayleigh | A=3, γ=2.5 |
| Tham số MGCC — Rician | A=5, γ=1.5 |
| Cấu hình anten | SISO, SIMO 1×2, MIMO 2×4, 2×8 |
| Dải SNR | 0–30 dB (bước 5 dB) |
| Số khung Monte Carlo | 500 khung/điểm SNR |

---

## Kết quả chính

| Kịch bản | Kết quả |
|----------|---------|
| CCDF PAPR | Giảm ~4–5 dB tại CCDF = 10⁻³ |
| MSE SISO, L=16 | Proposed-Transmitter ≈ baseline (Np = L, không có dư thừa ràng buộc) |
| MSE SISO, L=8,4 | Proposed-Transmitter-Rician vượt trội rõ rệt (MSE ~0.02 vs ~0.15) |
| BER SIMO 1×2 | Proposed-Transmitter-Rician đạt ~4×10⁻⁵ ở 30 dB |
| BER MIMO 2×4, 2×8 | Bão hòa cao do thiếu pilot trực giao giữa các anten phát |
| Tương quan anten ρ | Hiệu năng suy giảm rõ khi ρ > 0.6, mất phân tập khi ρ = 1 |

---

## Lưu ý kỹ thuật

**Về công thức MGCC:**  
Hàm `mgcc_expand.m` dùng nghịch đảo đại số đúng của hàm nén:
```
|x| = ( |y| / (A^γ - |y|) )^(1/γ)
```
khác với công thức (6) trong một số tài liệu tham khảo (đã kiểm chứng bằng giải tích, công thức gốc cho kết quả sai khi |y| > A).

**Về cấu hình MIMO đa anten phát:**  
Kết quả MIMO 2×4 và 2×8 bị bão hòa BER cao do mô hình mô phỏng hiện tại chưa sử dụng pilot trực giao giữa các anten phát. Đây là hướng cải tiến cho nghiên cứu tiếp theo.

---

## Tài liệu tham khảo

1. I. Khan, M. M. Hasan, and M. Cheffena, "Transmitter-assisted joint data-aided channel estimation and PAPR reduction scheme in wireless fading channels," *Scientific Reports*, vol. 16, no. 8015, 2026.
2. M. M. Hasan and M. M. H. Foad, "Modified gamma correction companding for PAPR reduction in OFDM systems," *Circuits, Systems, and Signal Processing*, vol. 37, pp. 4431–4454, 2018.
3. I. Khan, M. Cheffena, and M. M. Hasan, "Data aided channel estimation for MIMO-OFDM wireless systems using reliable carriers," *IEEE Access*, vol. 11, pp. 47836–47847, 2023.

---
