clear; close all; clc;

% =========================================================================
%  Fig 4 - CCDF của PAPR: General OFDM vs MGCC-OFDM
%  (Rayleigh: A=3, gamma=2.5  |  Rician: A=5, gamma=1.5)
%
%  File này KHÔNG gọi Test()/Test_1() vì PAPR là đặc tính của tín hiệu
%  phát (miền thời gian sau IFFT/MGCC), không phụ thuộc kênh hay SNR.
%  Dùng lại đúng tham số N_fft, Np, pilot_sym như trong Test.m để đảm
%  bảo nhất quán.
% =========================================================================

fprintf('\n=== Fig 4: CCDF cua PAPR ===\n');

%% ----------------------------------------------------------------
%  Tham số hệ thống (khớp với Test.m)
%% ----------------------------------------------------------------
N_fft   = 256;
Np      = 16;
Nd      = N_fft - Np;
M       = 4;                       % 4QAM cho khảo sát PAPR (đại diện)
num_frames = 5000;                 % số khung OFDM độc lập để vẽ CCDF mượt

pilot_interval = floor(N_fft / Np);
pilot_index    = 1 : pilot_interval : N_fft;
data_index     = setdiff(1:N_fft, pilot_index);
pilot_sym      = (1 + 1j) / sqrt(2);

% Tham số MGCC theo từng kênh (đúng như trong Test.m)
A_ray = 3;   gamma_ray = 2.5;   % Rayleigh
A_ric = 5;   gamma_ric = 1.5;   % Rician

%% ----------------------------------------------------------------
%  Tính PAPR cho từng trường hợp: General OFDM, MGCC-Rayleigh, MGCC-Rician
%% ----------------------------------------------------------------
papr_general = zeros(1, num_frames);
papr_mgcc_ray = zeros(1, num_frames);
papr_mgcc_ric = zeros(1, num_frames);

for fr = 1:num_frames
    data_bits = randi([0 1], log2(M)*Nd, 1);
    mod_syms  = qammod(data_bits, M, 'UnitAveragePower', true, 'InputType', 'bit');

    X = zeros(1, N_fft);
    X(pilot_index) = pilot_sym;
    X(data_index)  = mod_syms;

    x_time = ifft(X, N_fft) * sqrt(N_fft);

    % --- General OFDM (không companding) ---
    p_avg_gen  = mean(abs(x_time).^2);
    p_peak_gen = max(abs(x_time).^2);
    papr_general(fr) = p_peak_gen / p_avg_gen;

    % --- MGCC-OFDM, tham số Rayleigh ---
    x_mgcc_ray = mgcc_compress(x_time, A_ray, gamma_ray);
    p_scale_ray = sqrt(mean(abs(x_time(:)).^2) / mean(abs(x_mgcc_ray(:)).^2));
    x_mgcc_ray = x_mgcc_ray * p_scale_ray;
    p_avg  = mean(abs(x_mgcc_ray).^2);
    p_peak = max(abs(x_mgcc_ray).^2);
    papr_mgcc_ray(fr) = p_peak / p_avg;

    % --- MGCC-OFDM, tham số Rician ---
    x_mgcc_ric = mgcc_compress(x_time, A_ric, gamma_ric);
    p_scale_ric = sqrt(mean(abs(x_time(:)).^2) / mean(abs(x_mgcc_ric(:)).^2));
    x_mgcc_ric = x_mgcc_ric * p_scale_ric;
    p_avg  = mean(abs(x_mgcc_ric).^2);
    p_peak = max(abs(x_mgcc_ric).^2);
    papr_mgcc_ric(fr) = p_peak / p_avg;
end

%% ----------------------------------------------------------------
%  Tính CCDF:  P(PAPR > PAPR0)
%% ----------------------------------------------------------------
papr0_dB = 0:0.25:12;             % trục hoành (dB)
papr0    = 10.^(papr0_dB/10);

ccdf_general  = zeros(size(papr0));
ccdf_mgcc_ray = zeros(size(papr0));
ccdf_mgcc_ric = zeros(size(papr0));

for k = 1:length(papr0)
    ccdf_general(k)  = mean(papr_general  > papr0(k));
    ccdf_mgcc_ray(k) = mean(papr_mgcc_ray > papr0(k));
    ccdf_mgcc_ric(k) = mean(papr_mgcc_ric > papr0(k));
end

%% ----------------------------------------------------------------
%  Vẽ Fig 4
%% ----------------------------------------------------------------
figure('Name','Fig4'); hold on; grid on;
semilogy(papr0_dB, ccdf_general,  'k-',  'LineWidth',1.8, 'DisplayName','General OFDM');
semilogy(papr0_dB, ccdf_mgcc_ray, 'b--', 'LineWidth',1.8, 'DisplayName','MGCC-OFDM (Rayleigh: A=3, \\gamma=2.5)');
semilogy(papr0_dB, ccdf_mgcc_ric, 'r-.', 'LineWidth',1.8, 'DisplayName','MGCC-OFDM (Rician: A=5, \\gamma=1.5)');
xlabel('PAPR_0 (dB)'); ylabel('CCDF: P(PAPR > PAPR_0)');
legend('Location','southwest');
set(gca,'YScale','log');
ylim([1e-3 1]);
xlim([0 12]);
grid on; grid minor;
set(gca,'FontSize',11);
set(gca,'FontName','Times New Roman');
box on;
title('Fig 4 - CCDF cua PAPR: General OFDM vs MGCC-OFDM');

%% ----------------------------------------------------------------
%  In ra mức giảm PAPR tại CCDF = 1e-3 (mốc tham chiếu thường dùng)
%% ----------------------------------------------------------------
target_ccdf = 1e-3;
[~, idx_gen] = min(abs(ccdf_general  - target_ccdf));
[~, idx_ray] = min(abs(ccdf_mgcc_ray - target_ccdf));
[~, idx_ric] = min(abs(ccdf_mgcc_ric - target_ccdf));

fprintf('\nTai CCDF = 1e-3:\n');
fprintf('  General OFDM     : PAPR0 = %.2f dB\n', papr0_dB(idx_gen));
fprintf('  MGCC-OFDM Rayleigh: PAPR0 = %.2f dB  (giam %.2f dB)\n', ...
        papr0_dB(idx_ray), papr0_dB(idx_gen)-papr0_dB(idx_ray));
fprintf('  MGCC-OFDM Rician   : PAPR0 = %.2f dB  (giam %.2f dB)\n', ...
        papr0_dB(idx_ric), papr0_dB(idx_gen)-papr0_dB(idx_ric));

% Lưu số liệu để dùng lại (vẽ bằng công cụ khác nếu cần, hoặc đối chiếu)
save('-text', 'fig4_papr_data.txt', 'papr0_dB', 'ccdf_general', 'ccdf_mgcc_ray', 'ccdf_mgcc_ric');
fprintf('\nDa luu so lieu vao fig4_papr_data.txt\n');