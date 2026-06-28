clear; close all; clc;

SNR_dB = 0:5:30;
N      = length(SNR_dB);

% =========================================================================
%  FILE NÀY BỔ SUNG Fig 6-15 CHO Plot_test_1.m (vốn chỉ có Fig 5)
%  Toàn bộ dùng Test_1.m (công thức MGCC dạng lũy thừa căn, KHÁC với
%  Test.m). Cấu trúc hình giống hệt Plot_test.m để dễ đối chiếu.
%  Chạy Plot_test_1.m (Fig 5) riêng trước, sau đó chạy file này.
% =========================================================================


%% =========================================================================
%  Fig 5 – MSE, SISO, LS estimator, BPSK (M=2), L=16
%% =========================================================================
fprintf('\n=== Fig 5: SISO LS BPSK L=16 ===\n');

mse_base_ray_ls5 = zeros(1,N);  mse_tx_ray_ls5 = zeros(1,N);  mse_rx_ray_ls5 = zeros(1,N);
mse_base_ric_ls5 = zeros(1,N);  mse_tx_ric_ls5 = zeros(1,N);  mse_rx_ric_ls5 = zeros(1,N);

for i = 1:N
    [mse_base_ray_ls5(i), mse_tx_ray_ls5(i), mse_rx_ray_ls5(i), ~] = ...
        Test(SNR_dB(i), 'Rayleigh', 'LS', 16, 2, 0);
    [mse_base_ric_ls5(i), mse_tx_ric_ls5(i), mse_rx_ric_ls5(i), ~] = ...
        Test(SNR_dB(i), 'Rician',   'LS', 16, 2, 0);
end

figure('Name','Fig5'); hold on; grid on;
semilogy(SNR_dB, mse_base_ray_ls5, 'g-',   'LineWidth',1.8, 'DisplayName','LS-Rayleigh');
semilogy(SNR_dB, mse_tx_ray_ls5,   'b--',  'LineWidth',1.8, 'DisplayName','LS-Proposed-Transmitter-Rayleigh');
semilogy(SNR_dB, mse_tx_ric_ls5,   'b-p',  'LineWidth',1.8, 'DisplayName','LS-Proposed-Transmitter-Rician');
semilogy(SNR_dB, mse_base_ric_ls5, 'g-o',  'LineWidth',1.8, 'DisplayName','LS-Rician');
semilogy(SNR_dB, mse_rx_ray_ls5,   'm--s', 'LineWidth',1.8, 'DisplayName','LS-Proposed-Receiver-Rayleigh');
semilogy(SNR_dB, mse_rx_ric_ls5,   'm-*',  'LineWidth',1.8, 'DisplayName','LS-Proposed-Receiver-Rician');
xlabel('SNR (dB)'); ylabel('MSE'); legend('Location','southwest');
set(gca,'YScale','log');
ylim([1e-3 1e2]);
xlim([0 30]);
xticks(0:5:30);
grid on;
grid minor;
set(gca,'FontSize',11);
set(gca,'FontName','Times New Roman');
box on;
title('Fig 5 – MSE: SISO LS BPSK L=16');




%% =========================================================================
%  Fig 6 – MSE, SISO, LMMSE estimator, BPSK (M=2), L=16   (Test_1)
%% =========================================================================
fprintf('\n=== Fig 6 (Test_1): SISO LMMSE BPSK L=16 ===\n');

mse_base_ray_mm6 = zeros(1,N);  mse_tx_ray_mm6 = zeros(1,N);  mse_rx_ray_mm6 = zeros(1,N);
mse_base_ric_mm6 = zeros(1,N);  mse_tx_ric_mm6 = zeros(1,N);  mse_rx_ric_mm6 = zeros(1,N);

for i = 1:N
    [mse_base_ray_mm6(i), mse_tx_ray_mm6(i), mse_rx_ray_mm6(i), ~] = ...
        Test_1(SNR_dB(i), 'Rayleigh', 'LMMSE', 16, 2, 0);
    [mse_base_ric_mm6(i), mse_tx_ric_mm6(i), mse_rx_ric_mm6(i), ~] = ...
        Test_1(SNR_dB(i), 'Rician',   'LMMSE', 16, 2, 0);
end

figure('Name','Fig6_Test1'); hold on; grid on;
semilogy(SNR_dB, mse_base_ray_mm6, 'g-',   'LineWidth',1.8, 'DisplayName','LMMSE-Rayleigh');
semilogy(SNR_dB, mse_tx_ray_mm6,   'b--',  'LineWidth',1.8, 'DisplayName','LMMSE-Proposed-Transmitter-Rayleigh');
semilogy(SNR_dB, mse_tx_ric_mm6,   'b-p',  'LineWidth',1.8, 'DisplayName','LMMSE-Proposed-Transmitter-Rician');
semilogy(SNR_dB, mse_base_ric_mm6, 'g-o',  'LineWidth',1.8, 'DisplayName','LMMSE-Rician');
semilogy(SNR_dB, mse_rx_ray_mm6,   'm--s', 'LineWidth',1.8, 'DisplayName','LMMSE-Proposed-Receiver-Rayleigh');
semilogy(SNR_dB, mse_rx_ric_mm6,   'm-*',  'LineWidth',1.8, 'DisplayName','LMMSE-Proposed-Receiver-Rician');
xlabel('SNR (dB)'); ylabel('MSE'); legend('Location','southwest');
set(gca,'YScale','log'); ylim([1e-3 1e2]); xlim([0 30]); xticks(0:5:30);
grid on; grid minor; set(gca,'FontSize',11); set(gca,'FontName','Times New Roman');
box on; title('Fig 6 (Test_1) – MSE: SISO LMMSE BPSK L=16');

%% =========================================================================
%  Fig 7 – MSE, SISO, TX-assisted DACE, 4QAM (M=4), L=8   (Test_1)
%% =========================================================================
fprintf('\n=== Fig 7 (Test_1): SISO TX-DACE 4QAM L=8 ===\n');

mse_b7_ls  = zeros(1,N); mse_tx7_ray_ls  = zeros(1,N); mse_tx7_ric_ls  = zeros(1,N);
mse_b7_mm  = zeros(1,N); mse_tx7_ray_mm  = zeros(1,N); mse_tx7_ric_mm  = zeros(1,N);

for i = 1:N
    [mse_b7_ls(i),  mse_tx7_ray_ls(i),  ~, ~] = Test_1(SNR_dB(i),'Rayleigh','LS',   8, 4, 0);
    [~,             mse_tx7_ric_ls(i),   ~, ~] = Test_1(SNR_dB(i),'Rician',  'LS',   8, 4, 0);
    [mse_b7_mm(i),  mse_tx7_ray_mm(i),  ~, ~] = Test_1(SNR_dB(i),'Rayleigh','LMMSE',8, 4, 0);
    [~,             mse_tx7_ric_mm(i),   ~, ~] = Test_1(SNR_dB(i),'Rician',  'LMMSE',8, 4, 0);
end

figure('Name','Fig7_Test1'); hold on; grid on;
semilogy(SNR_dB, mse_b7_ls,        'b--',  'LineWidth',1.8, 'DisplayName','LS');
semilogy(SNR_dB, mse_tx7_ray_ls,   'g-',   'LineWidth',1.8, 'DisplayName','LS-Proposed-Transmitter-Rayleigh');
semilogy(SNR_dB, mse_tx7_ric_ls,   'm-p',  'LineWidth',1.8, 'DisplayName','LS-Proposed-Transmitter-Rician');
semilogy(SNR_dB, mse_b7_mm,        'b--s', 'LineWidth',1.8, 'DisplayName','LMMSE');
semilogy(SNR_dB, mse_tx7_ray_mm,   'k--',  'LineWidth',1.8, 'DisplayName','LMMSE-Proposed-Transmitter-Rayleigh');
semilogy(SNR_dB, mse_tx7_ric_mm,   'm-*',  'LineWidth',1.8, 'DisplayName','LMMSE-Proposed-Transmitter-Rician');
xlabel('SNR (dB)'); ylabel('MSE'); legend('Location','southwest');
set(gca,'YScale','log'); ylim([1e-3 1e2]); xlim([0 30]); xticks(0:5:30);
grid on; grid minor; set(gca,'FontSize',11); set(gca,'FontName','Times New Roman');
box on; title('Fig 7 (Test_1) – MSE: SISO TX-DACE 4QAM L=8');

%% =========================================================================
%  Fig 8 – MSE, SISO, TX-assisted DACE, 8PSK (M=8), L=8   (Test_1)
%% =========================================================================
fprintf('\n=== Fig 8 (Test_1): SISO TX-DACE 8PSK L=8 ===\n');

mse_b8_ls  = zeros(1,N); mse_tx8_ray_ls  = zeros(1,N); mse_tx8_ric_ls  = zeros(1,N);
mse_b8_mm  = zeros(1,N); mse_tx8_ray_mm  = zeros(1,N); mse_tx8_ric_mm  = zeros(1,N);

for i = 1:N
    [mse_b8_ls(i),  mse_tx8_ray_ls(i),  ~, ~] = Test_1(SNR_dB(i),'Rayleigh','LS',   8, 8, 0);
    [~,             mse_tx8_ric_ls(i),   ~, ~] = Test_1(SNR_dB(i),'Rician',  'LS',   8, 8, 0);
    [mse_b8_mm(i),  mse_tx8_ray_mm(i),  ~, ~] = Test_1(SNR_dB(i),'Rayleigh','LMMSE',8, 8, 0);
    [~,             mse_tx8_ric_mm(i),   ~, ~] = Test_1(SNR_dB(i),'Rician',  'LMMSE',8, 8, 0);
end

figure('Name','Fig8_Test1'); hold on; grid on;
semilogy(SNR_dB, mse_b8_ls,        'b--',  'LineWidth',1.8, 'DisplayName','LS');
semilogy(SNR_dB, mse_tx8_ray_ls,   'g-',   'LineWidth',1.8, 'DisplayName','LS-Proposed-Transmitter-Rayleigh');
semilogy(SNR_dB, mse_tx8_ric_ls,   'm-p',  'LineWidth',1.8, 'DisplayName','LS-Proposed-Transmitter-Rician');
semilogy(SNR_dB, mse_b8_mm,        'b--s', 'LineWidth',1.8, 'DisplayName','LMMSE');
semilogy(SNR_dB, mse_tx8_ray_mm,   'k--',  'LineWidth',1.8, 'DisplayName','LMMSE-Proposed-Transmitter-Rayleigh');
semilogy(SNR_dB, mse_tx8_ric_mm,   'm-*',  'LineWidth',1.8, 'DisplayName','LMMSE-Proposed-Transmitter-Rician');
xlabel('SNR (dB)'); ylabel('MSE'); legend('Location','southwest');
set(gca,'YScale','log'); ylim([1e-3 1e2]); xlim([0 30]); xticks(0:5:30);
grid on; grid minor; set(gca,'FontSize',11); set(gca,'FontName','Times New Roman');
box on; title('Fig 8 (Test_1) – MSE: SISO TX-DACE 8PSK L=8');

%% =========================================================================
%  Fig 9 – MSE, SISO, TX-assisted DACE, 4QAM, L=4   (Test_1)
%% =========================================================================
fprintf('\n=== Fig 9 (Test_1): SISO TX-DACE 4QAM L=4 ===\n');

mse_b9_ls  = zeros(1,N); mse_tx9_ray_ls  = zeros(1,N); mse_tx9_ric_ls  = zeros(1,N);
mse_b9_mm  = zeros(1,N); mse_tx9_ray_mm  = zeros(1,N); mse_tx9_ric_mm  = zeros(1,N);

for i = 1:N
    [mse_b9_ls(i),  mse_tx9_ray_ls(i),  ~, ~] = Test_1(SNR_dB(i),'Rayleigh','LS',   4, 4, 0);
    [~,             mse_tx9_ric_ls(i),   ~, ~] = Test_1(SNR_dB(i),'Rician',  'LS',   4, 4, 0);
    [mse_b9_mm(i),  mse_tx9_ray_mm(i),  ~, ~] = Test_1(SNR_dB(i),'Rayleigh','LMMSE',4, 4, 0);
    [~,             mse_tx9_ric_mm(i),   ~, ~] = Test_1(SNR_dB(i),'Rician',  'LMMSE',4, 4, 0);
end

figure('Name','Fig9_Test1'); hold on; grid on;
semilogy(SNR_dB, mse_b9_ls,        'b--',  'LineWidth',1.8, 'DisplayName','LS');
semilogy(SNR_dB, mse_tx9_ray_ls,   'g-',   'LineWidth',1.8, 'DisplayName','LS-Proposed-Transmitter-Rayleigh');
semilogy(SNR_dB, mse_tx9_ric_ls,   'm-p',  'LineWidth',1.8, 'DisplayName','LS-Proposed-Transmitter-Rician');
semilogy(SNR_dB, mse_b9_mm,        'b--s', 'LineWidth',1.8, 'DisplayName','LMMSE');
semilogy(SNR_dB, mse_tx9_ray_mm,   'k--',  'LineWidth',1.8, 'DisplayName','LMMSE-Proposed-Transmitter-Rayleigh');
semilogy(SNR_dB, mse_tx9_ric_mm,   'm-*',  'LineWidth',1.8, 'DisplayName','LMMSE-Proposed-Transmitter-Rician');
xlabel('SNR (dB)'); ylabel('MSE'); legend('Location','southwest');
set(gca,'YScale','log'); ylim([1e-3 1e2]); xlim([0 30]); xticks(0:5:30);
grid on; grid minor; set(gca,'FontSize',11); set(gca,'FontName','Times New Roman');
box on; title('Fig 9 (Test_1) – MSE: SISO TX-DACE 4QAM L=4');

%% =========================================================================
%  Fig 10 – MSE, SISO, TX-assisted DACE, 8PSK, L=4   (Test_1)
%% =========================================================================
fprintf('\n=== Fig 10 (Test_1): SISO TX-DACE 8PSK L=4 ===\n');

mse_b10_ls  = zeros(1,N); mse_tx10_ray_ls = zeros(1,N); mse_tx10_ric_ls = zeros(1,N);
mse_b10_mm  = zeros(1,N); mse_tx10_ray_mm = zeros(1,N); mse_tx10_ric_mm = zeros(1,N);

for i = 1:N
    [mse_b10_ls(i),  mse_tx10_ray_ls(i),  ~, ~] = Test_1(SNR_dB(i),'Rayleigh','LS',   4, 8, 0);
    [~,              mse_tx10_ric_ls(i),   ~, ~] = Test_1(SNR_dB(i),'Rician',  'LS',   4, 8, 0);
    [mse_b10_mm(i),  mse_tx10_ray_mm(i),  ~, ~] = Test_1(SNR_dB(i),'Rayleigh','LMMSE',4, 8, 0);
    [~,              mse_tx10_ric_mm(i),   ~, ~] = Test_1(SNR_dB(i),'Rician',  'LMMSE',4, 8, 0);
end

figure('Name','Fig10_Test1'); hold on; grid on;
semilogy(SNR_dB, mse_b10_ls,        'b--',  'LineWidth',1.8, 'DisplayName','LS');
semilogy(SNR_dB, mse_tx10_ray_ls,   'g-',   'LineWidth',1.8, 'DisplayName','LS-Proposed-Transmitter-Rayleigh');
semilogy(SNR_dB, mse_tx10_ric_ls,   'm-p',  'LineWidth',1.8, 'DisplayName','LS-Proposed-Transmitter-Rician');
semilogy(SNR_dB, mse_b10_mm,        'b--s', 'LineWidth',1.8, 'DisplayName','LMMSE');
semilogy(SNR_dB, mse_tx10_ray_mm,   'k--',  'LineWidth',1.8, 'DisplayName','LMMSE-Proposed-Transmitter-Rayleigh');
semilogy(SNR_dB, mse_tx10_ric_mm,   'm-*',  'LineWidth',1.8, 'DisplayName','LMMSE-Proposed-Transmitter-Rician');
xlabel('SNR (dB)'); ylabel('MSE'); legend('Location','southwest');
set(gca,'YScale','log'); ylim([1e-3 1e2]); xlim([0 30]); xticks(0:5:30);
grid on; grid minor; set(gca,'FontSize',11); set(gca,'FontName','Times New Roman');
box on; title('Fig 10 (Test_1) – MSE: SISO TX-DACE 8PSK L=4');

%% =========================================================================
%  Fig 11 – BER, SIMO (Nt=1, Nr=2), 4QAM & 8PSK, L=16   (Test_1)
%% =========================================================================
fprintf('\n=== Fig 11 (Test_1): SIMO 1x2 BER L=16 ===\n');

ber_base_4qam  = zeros(1,N); ber_tx_ray_4qam = zeros(1,N); ber_tx_ric_4qam = zeros(1,N);
ber_base_8psk  = zeros(1,N); ber_tx_ray_8psk = zeros(1,N); ber_tx_ric_8psk = zeros(1,N);

for i = 1:N
    [~, ~, ~, ber_base_4qam(i)]  = Test_1(SNR_dB(i),'Rayleigh','LS', 16, 4, 0, 1, 2);
    [~, ~, ~, ber_tx_ray_4qam(i)]= Test_1(SNR_dB(i),'Rayleigh','LS', 16, 4, 0, 1, 2);
    [~, ~, ~, ber_tx_ric_4qam(i)]= Test_1(SNR_dB(i),'Rician',  'LS', 16, 4, 0, 1, 2);
    [~, ~, ~, ber_base_8psk(i)]  = Test_1(SNR_dB(i),'Rayleigh','LS', 16, 8, 0, 1, 2);
    [~, ~, ~, ber_tx_ray_8psk(i)]= Test_1(SNR_dB(i),'Rayleigh','LS', 16, 8, 0, 1, 2);
    [~, ~, ~, ber_tx_ric_8psk(i)]= Test_1(SNR_dB(i),'Rician',  'LS', 16, 8, 0, 1, 2);
end

figure('Name','Fig11_Test1'); hold on; grid on;
semilogy(SNR_dB, ber_base_8psk,    'b--s', 'LineWidth',1.8, 'DisplayName','8PSK');
semilogy(SNR_dB, ber_tx_ray_8psk,  'g-',   'LineWidth',1.8, 'DisplayName','8PSK-Proposed-Transmitter-Rayleigh');
semilogy(SNR_dB, ber_tx_ric_8psk,  'm-p',  'LineWidth',1.8, 'DisplayName','8PSK-Proposed-Transmitter-Rician');
semilogy(SNR_dB, ber_base_4qam,    'b--^', 'LineWidth',1.8, 'DisplayName','4QAM');
semilogy(SNR_dB, ber_tx_ray_4qam,  'k-',   'LineWidth',1.8, 'DisplayName','4QAM-Proposed-Transmitter-Rayleigh');
semilogy(SNR_dB, ber_tx_ric_4qam,  'm-*',  'LineWidth',1.8, 'DisplayName','4QAM-Proposed-Transmitter-Rician');
xlabel('SNR (dB)'); ylabel('BER'); legend('Location','southwest');
set(gca,'YScale','log'); ylim([1e-7 1e0]); xlim([0 30]); xticks(0:5:30);
grid on; grid minor; set(gca,'FontSize',11); set(gca,'FontName','Times New Roman');
title('Fig 11 (Test_1) – BER: SIMO 1×2 4QAM & 8PSK L=16');

%% =========================================================================
%  Fig 13 – BER, 2×4 MIMO, 4QAM & 8PSK, L=16   (Test_1)
%% =========================================================================
fprintf('\n=== Fig 13 (Test_1): 2x4 MIMO BER L=16 ===\n');

ber_b13_8psk     = zeros(1,N); ber_tx13_ray_8psk = zeros(1,N); ber_tx13_ric_8psk = zeros(1,N);
ber_b13_4qam     = zeros(1,N); ber_tx13_ray_4qam = zeros(1,N); ber_tx13_ric_4qam = zeros(1,N);

for i = 1:N
    [~,~,~, ber_b13_8psk(i)]      = Test_1(SNR_dB(i),'Rayleigh','LS',16, 8,0, 2,4);
    [~,~,~, ber_tx13_ray_8psk(i)] = Test_1(SNR_dB(i),'Rayleigh','LS',16, 8,0, 2,4);
    [~,~,~, ber_tx13_ric_8psk(i)] = Test_1(SNR_dB(i),'Rician',  'LS',16, 8,0, 2,4);
    [~,~,~, ber_b13_4qam(i)]      = Test_1(SNR_dB(i),'Rayleigh','LS',16, 4,0, 2,4);
    [~,~,~, ber_tx13_ray_4qam(i)] = Test_1(SNR_dB(i),'Rayleigh','LS',16, 4,0, 2,4);
    [~,~,~, ber_tx13_ric_4qam(i)] = Test_1(SNR_dB(i),'Rician',  'LS',16, 4,0, 2,4);
end

figure('Name','Fig13_Test1'); hold on; grid on;
semilogy(SNR_dB, ber_b13_8psk,      'b--s', 'LineWidth',1.8, 'DisplayName','8PSK');
semilogy(SNR_dB, ber_tx13_ray_8psk, 'g-',   'LineWidth',1.8, 'DisplayName','8PSK-Proposed-Transmitter-Rayleigh');
semilogy(SNR_dB, ber_tx13_ric_8psk, 'm-p',  'LineWidth',1.8, 'DisplayName','8PSK-Proposed-Transmitter-Rician');
semilogy(SNR_dB, ber_b13_4qam,      'b--^', 'LineWidth',1.8, 'DisplayName','4QAM');
semilogy(SNR_dB, ber_tx13_ray_4qam, 'k-',   'LineWidth',1.8, 'DisplayName','4QAM-Proposed-Transmitter-Rayleigh');
semilogy(SNR_dB, ber_tx13_ric_4qam, 'm-*',  'LineWidth',1.8, 'DisplayName','4QAM-Proposed-Transmitter-Rician');
xlabel('SNR (dB)'); ylabel('BER'); legend('Location','southwest');
set(gca,'YScale','log'); ylim([1e-3 1e0]); xlim([0 30]); xticks(0:5:30);
grid on; grid minor; set(gca,'FontSize',11); set(gca,'FontName','Times New Roman');
title('Fig 13 (Test_1) – BER: 2×4 MIMO 4QAM & 8PSK L=16');

%% =========================================================================
%  Fig 14 – BER, 2×8 MIMO, 4QAM & 8PSK, L=16   (Test_1)
%% =========================================================================
fprintf('\n=== Fig 14 (Test_1): 2x8 MIMO BER L=16 ===\n');

ber_b14_8psk     = zeros(1,N); ber_tx14_ray_8psk = zeros(1,N); ber_tx14_ric_8psk = zeros(1,N);
ber_b14_4qam     = zeros(1,N); ber_tx14_ray_4qam = zeros(1,N); ber_tx14_ric_4qam = zeros(1,N);

for i = 1:N
    [~,~,~, ber_b14_8psk(i)]      = Test_1(SNR_dB(i),'Rayleigh','LS',16, 8,0, 2,8);
    [~,~,~, ber_tx14_ray_8psk(i)] = Test_1(SNR_dB(i),'Rayleigh','LS',16, 8,0, 2,8);
    [~,~,~, ber_tx14_ric_8psk(i)] = Test_1(SNR_dB(i),'Rician',  'LS',16, 8,0, 2,8);
    [~,~,~, ber_b14_4qam(i)]      = Test_1(SNR_dB(i),'Rayleigh','LS',16, 4,0, 2,8);
    [~,~,~, ber_tx14_ray_4qam(i)] = Test_1(SNR_dB(i),'Rayleigh','LS',16, 4,0, 2,8);
    [~,~,~, ber_tx14_ric_4qam(i)] = Test_1(SNR_dB(i),'Rician',  'LS',16, 4,0, 2,8);
end

figure('Name','Fig14_Test1'); hold on; grid on;
semilogy(SNR_dB, ber_b14_8psk,      'b--s', 'LineWidth',1.8, 'DisplayName','8PSK');
semilogy(SNR_dB, ber_tx14_ray_8psk, 'g-',   'LineWidth',1.8, 'DisplayName','8PSK-Proposed-Transmitter-Rayleigh');
semilogy(SNR_dB, ber_tx14_ric_8psk, 'm-p',  'LineWidth',1.8, 'DisplayName','8PSK-Proposed-Transmitter-Rician');
semilogy(SNR_dB, ber_b14_4qam,      'b--^', 'LineWidth',1.8, 'DisplayName','4QAM');
semilogy(SNR_dB, ber_tx14_ray_4qam, 'k-',   'LineWidth',1.8, 'DisplayName','4QAM-Proposed-Transmitter-Rayleigh');
semilogy(SNR_dB, ber_tx14_ric_4qam, 'm-*',  'LineWidth',1.8, 'DisplayName','4QAM-Proposed-Transmitter-Rician');
xlabel('SNR (dB)'); ylabel('BER'); legend('Location','southwest');
set(gca,'YScale','log'); ylim([1e-3 1e0]); xlim([0 30]); xticks(0:5:30);
grid on; grid minor; set(gca,'FontSize',11); set(gca,'FontName','Times New Roman');
title('Fig 14 (Test_1) – BER: 2×8 MIMO 4QAM & 8PSK L=16');

%% =========================================================================
%  Fig 15 – BER vs. receive antenna correlation rho, SISO, 4QAM   (Test_1)
%  SNR range 8:2:30
%% =========================================================================
fprintf('\n=== Fig 15 (Test_1): Correlation effect rho ===\n');

SNR15 = 8:2:30;
rho_vals = [0, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0];
ber15 = zeros(length(rho_vals), length(SNR15));

for ri = 1:length(rho_vals)
    for si = 1:length(SNR15)
        [~,~,~, ber15(ri,si)] = Test_1(SNR15(si),'Rician','LS',16,4, rho_vals(ri), 1, 2);
    end
    fprintf('  rho=%.1f done\n', rho_vals(ri));
end

figure('Name','Fig15_Test1'); hold on; grid on;
colors15 = lines(length(rho_vals));
markers  = {'p-','o-','s-','d-','^-','v-','*-'};
for ri = 1:length(rho_vals)
    semilogy(SNR15, ber15(ri,:), markers{ri}, ...
             'Color', colors15(ri,:), 'LineWidth',1.5, ...
             'DisplayName', sprintf('\\rho = %.1f', rho_vals(ri)));
end
xlabel('SNR (dB)'); ylabel('BER'); legend('Location','southwest');
set(gca,'YScale','log'); ylim([1e-6 1e-1]); xlim([0 30]); xticks(8:2:30);
grid on; grid minor; set(gca,'FontSize',11); set(gca,'FontName','Times New Roman');
title('Fig 15 (Test_1) – BER under varying antenna correlation (\rho)');

fprintf('\n=== Da ve xong Fig 6-15 (Test_1) ===\n');