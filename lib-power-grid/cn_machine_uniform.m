function colored_noise=cn_machine_uniform(N,b,w,fs)
f = fs * (0:N/2-1) / N;    % frequency

white_noise = w*rand(1, N);  

Y = fft(white_noise);

% create the 1/f relationship in frequency domain
f=[f,f(end)+fs/2*N];
Y(2:N/2+1) = Y(2:N/2+1) .* (1 ./ (f(2:N/2+1).^(b/2))); 
Y(N/2+2:N) = flip(conj(Y(2:N/2))); 

% ifft
colored_noise = ifft(Y);
colored_noise = colored_noise-0.5*w;

end
