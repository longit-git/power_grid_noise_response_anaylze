function result = color_gradient(beginning, ending, num)
% color_gradient 在两种 RGB 颜色之间线性插值，生成 num 行渐变色
%
% 用法：
%   colors = color_gradient([0.90 0.80 0.40], [0.40 0.70 0.40], 10);

result = zeros(num, length(beginning));
for i = 1:num
    result(i, :) = beginning + ((i-1)/num) * (ending - beginning);
end
end
