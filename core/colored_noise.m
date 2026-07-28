function noise = colored_noise(N, b, w, Fs)
% COLORED_NOISE  生成 1/f^b 有色噪声（频域滤波法）
%
% 输入：
%   N    — 输出长度（采样点数）
%   b    — 噪声颜色指数（0=白噪声, 1=粉噪声, 2=布朗噪声）
%   w    — 均匀白噪声幅度
%   Fs   — 采样率 [Hz]
%
% 输出：
%   noise — N×1 实数列向量

f = Fs * (0 : N/2-1) / N;          % 单边频率轴（含零频；下方滤波从第 2 点开始，跳过零频）
white = w * rand(1, N);
Y = fft(white);

f_full = [f, f(end) + Fs/(2*N)];   % 补全至 N/2+1 个频点
Y(2 : N/2+1) = Y(2 : N/2+1) ./ (f_full(2 : N/2+1) .^ (b/2));
Y(N/2+2 : N) = flip(conj(Y(2 : N/2)));

noise = real(ifft(Y));
noise = noise - 0.5*w;              % 去均值
noise = noise(:);                   % 确保列向量
end
