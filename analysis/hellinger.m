function d = hellinger(P, Q)
% hellinger Hellinger distance，输入为 probability 向量（和为 1）
%
%   d = ||sqrt(P) - sqrt(Q)||_2 / sqrt(2)

d = sqrt(sum((sqrt(P) - sqrt(Q)).^2)) / sqrt(2);
end
