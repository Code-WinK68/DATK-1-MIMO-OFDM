function xe = mgcc_expand(y, A, gamma)
% MGCC_EXPAND  Giãn tín hiệu - nghịch đảo của mgcc_compress
%
%   xe = mgcc_expand(y, A, gamma)
%
%   LƯU Ý QUAN TRỌNG: Công thức (6) trong bài báo gốc viết
%       H_gamma^-1(y) = sgn(y) * (A/|y| - 1)^(-1/gamma)
%   KHÔNG PHẢI là nghịch đảo đại số đúng của công thức (5). Đã kiểm
%   chứng bằng giải tích (sympy.solve): nghịch đảo đúng của
%       y = (A|x|)^gamma / (1+|x|^gamma)
%   là
%       |x| = ( |y| / (A^gamma - |y|) )^(1/gamma),   0 <= |y| < A^gamma
%   Công thức (6) của bài báo cho giá trị PHỨC (sai) khi |y| > A. Đây
%   là lỗi/typo trong bài báo gốc. Hàm này dùng công thức đại số ĐÚNG,
%   đã kiểm chứng round-trip với mgcc_compress (sai số ~1e-16).
%
%   GHI CHÚ THỰC NGHIỆM: qua kênh fading + AWGN, một tỉ lệ nhỏ (~1-2%)
%   mẫu |y| có thể vượt asymptote A^gamma do hệ số khuếch đại kênh.
%   Đây là hạn chế vật lý cố hữu của companding qua kênh fading (không
%   phải lỗi code). Áp dụng soft-clamp margin=0.97 để tránh nổ số.
%
%   Inputs:
%     y     : tín hiệu phức cần giãn (đã qua kênh + nhiễu)
%     A     : tham số biên độ (phải khớp với A dùng lúc nén)
%     gamma : tham số mức nén  (phải khớp với gamma dùng lúc nén)
%
%   Output:
%     xe    : tín hiệu đã giãn (ước lượng của tín hiệu gốc trước khi nén)

    yabs = abs(y);
    s    = y ./ max(yabs, eps);
    Agam = A^gamma;
    margin = 0.97;
    yabs_clamped = min(yabs, Agam * margin);
    ratio = yabs_clamped ./ max(Agam - yabs_clamped, eps);
    xe = s .* ratio.^(1/gamma);
end