function xc = mgcc_compress(x, A, gamma)
% MGCC_COMPRESS  Nén tín hiệu theo công thức (5) của bài báo gốc
%                (Khan, Hasan, Cheffena - Scientific Reports 2026)
%
%   xc = mgcc_compress(x, A, gamma)
%
%   H_gamma(x) = sgn(x) * (A|x|)^gamma / (1 + |x|^gamma)
%
%   Inputs:
%     x     : tín hiệu phức cần nén (vector/ma trận bất kỳ)
%     A     : tham số biên độ (Rayleigh: A=3, Rician: A=5)
%     gamma : tham số mức nén  (Rayleigh: gamma=2.5, Rician: gamma=1.5)
%
%   Output:
%     xc    : tín hiệu đã nén, |xc| -> A^gamma khi |x| -> vô cùng

    xabs = abs(x);
    s    = x ./ max(xabs, eps);          % "dấu phức" của x, |s| = 1
    num  = (A .* xabs).^gamma;
    den  = 1 + xabs.^gamma;
    xc   = s .* (num ./ den);
end