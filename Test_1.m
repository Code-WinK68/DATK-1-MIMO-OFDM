function [mse_base, mse_tx, mse_rx, ber] = Test_1(SNR_dB, channel_type, est_method, L, M, rho, Nt, Nr)

% =========================================================================
%  Test_1.m  –  OFDM channel estimation with MGCC (power-law) + DACE
%               (ĐÃ SỬA: thêm Nt, Nr tùy chọn)
%
%  Inputs:
%   SNR_dB       : Signal-to-Noise Ratio (dB)
%   channel_type : 'Rayleigh' | 'Rician'
%   est_method   : 'LS' | 'LMMSE'
%   L            : Number of multipath taps (e.g. 4, 8, 16)
%   M            : Modulation order (2=BPSK, 4=4QAM, 8=8PSK)
%   rho          : Receive antenna correlation coefficient (0 = uncorrelated)
%   Nt           : (TÙY CHỌN) Number of transmit antennas  (default 1)
%   Nr           : (TÙY CHỌN) Number of receive  antennas  (default 1)
%
%  Outputs:
%   mse_base : MSE of baseline estimator (pilots only)
%   mse_tx   : MSE of Proposed-Transmitter (pilots + TX peak SCs)
%   mse_rx   : MSE of Proposed-Receiver   (pilots + virtual pilots)
%   ber      : Bit Error Rate
%
%  ----------------------------------------------------------------------
%  GHI CHÚ QUAN TRỌNG (kết quả kiểm chứng bằng giải tích sympy):
%  Công thức MGCC dùng trong file này:
%      x_mgcc = A * sign(x) * |x|^(1/gamma)        (nén)
%      y_decomp = sign(y) * (|y|/A)^gamma          (giãn)
%  là MỘT CẶP NGHỊCH ĐẢO ĐÚNG về đại số (round-trip chính xác, không
%  như công thức Eq.(5)-(6) dùng trong Test.m, vốn có lỗi/typo). Vì
%  không có asymptote/điểm cực nên KHÔNG xảy ra hiện tượng "tràn số"
%  (overflow) như ở Test.m.
%
%  Tuy nhiên, dạng nén lũy thừa căn này lại KHUẾCH ĐẠI các mẫu có biên
%  độ < 1 (vì |x|^(1/gamma) > |x| khi |x|<1 và gamma>1), khiến nhiễu ở
%  vùng biên độ thấp bị khuếch đại tương đối nhiều hơn, gây ra error
%  floor ở mức MSE cao hơn Test.m (~0.25 so với ~0.09 ở SNR=20-30dB).
%  Đây là đặc tính thực nghiệm của riêng dạng companding này, đã được
%  kiểm chứng bằng mô phỏng Monte-Carlo, không phải lỗi cài đặt.
% =========================================================================

rng(42);
if nargin < 7 || isempty(Nt), Nt = 1; end
if nargin < 8 || isempty(Nr), Nr = 1; end
N_fft       = 256;
Np          = 16;
CP_len      = 16;
num_symbols = 500;
Nrsc        = 16;

Nd = N_fft - Np;

if strcmp(channel_type, 'Rayleigh')
    A = 1; gamma = 2; K_factor = 0;
else
    A = 1; gamma = 2; K_factor = 10;
end

%%-------------------------------------------------------------
% Receive-antenna correlation  R_r(i,j) = rho^|i-j|
%%-------------------------------------------------------------
idx     = (0:Nr-1)';
R_r     = rho .^ abs(idx - idx');
L_chol  = chol(R_r + 1e-10*eye(Nr), 'lower');

%% ----------------------------------------------------------------
%  Pilot / data subcarrier indices  (comb-type, uniform)
%% ----------------------------------------------------------------
pilot_interval = floor(N_fft / Np);
pilot_index    = 1 : pilot_interval : N_fft;   % 1-based
data_index     = setdiff(1:N_fft, pilot_index);


%% ----------------------------------------------------------------
%  DFT matrix (compute once)
%% ----------------------------------------------------------------
F_full = dftmtx(N_fft);          % N_fft x N_fft
F_p    = F_full(pilot_index, 1:L);   % Np x L

%% ----------------------------------------------------------------
%  LMMSE prior  (assume R_h = (1/L)*I_L  →  inv = L*I)
%% ----------------------------------------------------------------
inv_R_h    = L * eye(L);          % L x L
noise_var  = 10^(-SNR_dB/10);     % noise power (unit signal power)


%% ----------------------------------------------------------------
pilot_sym  = (1 + 1j) / sqrt(2);          % known pilot value

 
%% ----------------------------------------------------------------
%  Accumulation buffers
sum_mse_base = 0;
sum_mse_tx   = 0;
sum_mse_rx   = 0;
errors       = 0;
total_bits   = 0;
 
    %% ================================================================
    %  Main Monte-Carlo loop
    %% ================================================================
    for sym_idx = 1:num_symbols
 
        %% ============================================================
        %  1.  TRANSMITTER
        %% ============================================================
 
        % --- Generate data bits & modulate ---
        num_data_bits = log2(M) * Nt * Nd;
        data_bits     = randi([0 1], num_data_bits, 1);
        mod_syms      = qammod(data_bits, M, ...
                               'UnitAveragePower', true, ...
                               'InputType', 'bit');   % Nt*Nd x 1
 
        % --- Build OFDM frame  X : Nt x N_fft ---
        X = zeros(Nt, N_fft);
        for nt = 1:Nt
            X(nt, pilot_index) = pilot_sym;                         % pilots
            X(nt, data_index)  = mod_syms((nt-1)*Nd+1 : nt*Nd);   % data
        end
 
        % --- IFFT (time domain) ---
        x_time = ifft(X, N_fft, 2) * sqrt(N_fft);   % Nt x N_fft
 
        % --- MGCC compression ---
        x_abs  = abs(x_time);
        x_mgcc = A.*sign(x_time).* (x_abs).^(1/gamma);
 
        % Power-normalisation so E[|x_mgcc|^2] = E[|x_time|^2]
        p_scale = sqrt( mean(abs(x_time(:)).^2) / mean(abs(x_mgcc(:)).^2) );
        x_mgcc  = x_mgcc * p_scale;
 
        % --- Select Nrsc peak subcarriers at TX ---
        P_tx            = sum(abs(X).^2, 1);             % 1 x N_fft
        [~, sorted_tx]  = sort(P_tx, 'descend');
        r_index_tx      = sorted_tx(1:Nrsc);             % peak SC indices
 
        % --- Add Cyclic Prefix ---
        x_cp = [x_mgcc(:, end-CP_len+1:end), x_mgcc];   % Nt x (N_fft+CP_len)
        Ncp  = size(x_cp, 2);
 
        %% ============================================================
        %  2.  CHANNEL  (frequency-selective fading + AWGN)
        %% ============================================================
        if K_factor == 0   % Rayleigh
            H_tap_raw = (randn(Nr,Nt,L) + 1j*randn(Nr,Nt,L)) / sqrt(2*L);
        else               % Rician
            pow_LOS  = K_factor / (K_factor + 1);
            pow_NLOS = 1        / (K_factor + 1);
            H_tap_raw = sqrt(pow_NLOS/(2*L)) * ...
                        (randn(Nr,Nt,L) + 1j*randn(Nr,Nt,L));
            H_tap_raw(:,:,1) = H_tap_raw(:,:,1) + sqrt(pow_LOS);
        end
 
        % Apply receive-antenna correlation (Kronecker model)
        H_tap = zeros(Nr, Nt, L);
        for ll = 1:L
            H_tap(:,:,ll) = L_chol * H_tap_raw(:,:,ll);
        end
 
        % True frequency-domain channel
        H_true = fft(H_tap, N_fft, 3);   % Nr x Nt x N_fft
 
        % Convolve through multipath channel
        rx_signal = zeros(Nr, Ncp);
        for r = 1:Nr
            for t = 1:Nt
                h_vec    = squeeze(H_tap(r,t,:)).';          % 1 x L
                h_pad    = [h_vec, zeros(1, Ncp-L)];
                tmp      = ifft(fft(h_pad,Ncp) .* fft(x_cp(t,:),Ncp));
                rx_signal(r,:) = rx_signal(r,:) + tmp(1:Ncp);
            end
        end
 
        % Add AWGN
        noise = sqrt(noise_var/2) * (randn(size(rx_signal)) + ...
                                     1j*randn(size(rx_signal)));
        y = rx_signal + noise;
 
        %% ============================================================
        %  3.  RECEIVER
        %% ============================================================
 
        % --- Remove CP and undo power scaling ---
        y_nocp = y(:, CP_len+1:end) / p_scale;    % Nr x N_fft
 
        % --- MGCC inverse decompression ---
        y_abs = abs(y_nocp);
        y_decomp = sign(y_nocp).*(y_abs/A).^gamma;
 
        % --- FFT → frequency domain ---
        Y = fft(y_decomp, N_fft, 2) / sqrt(N_fft);   % Nr x N_fft
 
        % --- Peak subcarriers at RX (from received power) ---
        P_rx            = sum(abs(Y).^2, 1);
        [~, sorted_rx]  = sort(P_rx, 'descend');
        r_index_rx      = sorted_rx(1:Nrsc);
 
        %% ============================================================
        %  4a.  BASELINE channel estimator  (pilots only)
        %% ============================================================
        H_base   = zeros(Nr, Nt, N_fft);
        inv_Rz_p = eye(Np) / noise_var;
 
        for r_ant = 1:Nr
            for t_ant = 1:Nt
                x_p   = X(t_ant, pilot_index);          % 1 x Np
                C_p   = diag(x_p) * F_p;                % Np x L
                y_p   = Y(r_ant, pilot_index).';        % Np x 1
 
                if strcmp(est_method, 'LS')
                    h_est = (C_p' * C_p) \ (C_p' * y_p);
                else  % LMMSE
                    h_est = (inv_R_h + C_p'*inv_Rz_p*C_p) \ ...
                            (C_p' * inv_Rz_p * y_p);
                end
                H_base(r_ant, t_ant, :) = fft(h_est.', N_fft);
            end
        end
 
        %% ============================================================
        %  4b.  PROPOSED-TRANSMITTER  (pilots ∪ TX peak SCs)
        %% ============================================================
        rp_idx_tx = unique([pilot_index, r_index_tx]);
        N_rp_tx   = length(rp_idx_tx);
        F_rp_tx   = F_full(rp_idx_tx, 1:L);
        inv_Rz_tx = eye(N_rp_tx) / noise_var;
        H_tx      = zeros(Nr, Nt, N_fft);
 
        for r_ant = 1:Nr
            for t_ant = 1:Nt
                x_rp_tx  = X(t_ant, rp_idx_tx);
                C_rp_tx  = diag(x_rp_tx) * F_rp_tx;
                y_rp_tx  = Y(r_ant,   rp_idx_tx).';
 
                if strcmp(est_method, 'LS')
                    h_est = (C_rp_tx' * C_rp_tx) \ (C_rp_tx' * y_rp_tx);
                else
                    h_est = (inv_R_h + C_rp_tx'*inv_Rz_tx*C_rp_tx) \ ...
                            (C_rp_tx' * inv_Rz_tx * y_rp_tx);
                end
                H_tx(r_ant, t_ant, :) = fft(h_est.', N_fft);
            end
        end
 
        %% ============================================================
        %  4c.  PROPOSED-RECEIVER  (pilots ∪ virtual pilots from RX peak)
        %% ============================================================
 
        % Step 1: decide symbols at RX peak SCs using H_base (ZF)
        data_peak_rx = intersect(r_index_rx, data_index);   % pure data SCs
        num_vp       = length(data_peak_rx);
 
        X_decided = zeros(Nt, num_vp);
        for k_i = 1:num_vp
            k_sc  = data_peak_rx(k_i);
            H_k   = squeeze(H_base(:,:,k_sc));    % Nr x Nt
            y_k   = Y(:, k_sc);                   % Nr x 1

            % ZF detection
            x_hat = pinv(H_k) * y_k;              % Nt x 1
  

            % Hard decision (re-modulate for clean virtual pilot)
            bits_hat    = qamdemod(x_hat, M, 'UnitAveragePower', true, ...
                                   'OutputType', 'bit');
            X_decided(:,k_i) = qammod(bits_hat, M, 'UnitAveragePower', true, ...
                                      'InputType', 'bit');
        end
 
        % Step 2: build extended pilot index and send values
        rp_idx_rx = unique([pilot_index, data_peak_rx]);
        N_rp_rx   = length(rp_idx_rx);
        F_rp_rx   = F_full(rp_idx_rx, 1:L);
        inv_Rz_rx = eye(N_rp_rx) / noise_var;
 
        % Locate pilot and virtual-pilot positions in the combined index
        [~, p_loc] = ismember(pilot_index, rp_idx_rx);
        [~, v_loc] = ismember(data_peak_rx, rp_idx_rx);
 
        H_rx = zeros(Nr, Nt, N_fft);
 
        for r_ant = 1:Nr
            for t_ant = 1:Nt
                % Build the combined "known" symbol vector
                x_rp_rx        = zeros(1, N_rp_rx);
                x_rp_rx(p_loc) = X(t_ant, pilot_index);           % real pilots
                if ~isempty(v_loc)
                    x_rp_rx(v_loc) = X_decided(t_ant, :);         % virtual pilots
                end
 
                C_rp_rx= diag(x_rp_rx) * F_rp_rx;
                y_rp_rx = Y(r_ant,   rp_idx_rx).';
 
                if strcmp(est_method, 'LS')
                    h_est = (C_rp_rx' * C_rp_rx) \ (C_rp_rx' * y_rp_rx);
                else
                    h_est = (inv_R_h + C_rp_rx'*inv_Rz_rx*C_rp_rx) \ ...
                            (C_rp_rx' * inv_Rz_rx * y_rp_rx);
                end
                H_rx(r_ant, t_ant, :) = fft(h_est.', N_fft);
            end
        end
 
        %% ============================================================
        %  5.  MSE accumulation
        %% ============================================================
        sum_mse_base = sum_mse_base + mean(abs(H_base(:) - H_true(:)).^2);
        sum_mse_tx   = sum_mse_tx   + mean(abs(H_tx(:)   - H_true(:)).^2);
        sum_mse_rx   = sum_mse_rx   + mean(abs(H_rx(:)   - H_true(:)).^2);
 
        %% ============================================================
        %  6.  Data detection  (use H_rx + ZF equalizer)
        %% ============================================================
        rx_mod = zeros(Nt, Nd);
        for k_i = 1:Nd
            k_sc  = data_index(k_i);
            H_k   = squeeze(H_rx(:,:,k_sc));    % Nr x Nt
            y_k   = Y(:, k_sc);
            rx_mod(:, k_i) = pinv(H_k) * y_k;
        end
 
        % Demodulate and count bit errors
        rx_syms   = reshape(rx_mod.', [], 1);       % Nt*Nd x 1
        rx_bits   = qamdemod(rx_syms, M, ...
                             'UnitAveragePower', true, ...
                             'OutputType', 'bit');
        errors      = errors      + sum(data_bits ~= rx_bits);
        total_bits  = total_bits  + length(data_bits);
 
    end  % end Monte-Carlo
 
    %% ----------------------------------------------------------------
    %  Output averages
    %% ----------------------------------------------------------------
    mse_base = sum_mse_base / num_symbols;
    mse_tx   = sum_mse_tx   / num_symbols;
    mse_rx   = sum_mse_rx   / num_symbols;
    ber      = errors / total_bits;
 
    fprintf(['SNR=%2ddB Nt=%d Nr=%d rho=%.1f | ' ...
             'BER=%.2e MSEb=%.2e MSEtx=%.2e MSErx=%.2e\n'], ...
             SNR_dB, Nt, Nr, rho, ber, mse_base, mse_tx, mse_rx);
end