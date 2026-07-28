function p = gauss_pdf(x, mu, sigma)
% gauss_pdf Gaussian pdf N(mu, sigma^2) 在 x 处的取值
%
% 用法：
%   p = gauss_pdf(centers, mu, sigma)

p = (1 ./ (sigma * sqrt(2*pi))) .* exp(-0.5 * ((x - mu) ./ sigma).^2);
end
