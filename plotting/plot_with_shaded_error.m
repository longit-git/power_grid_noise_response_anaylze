function plot_with_shaded_error(x, y_mean, y_err, color)
% plot_with_shaded_error 绘制均值曲线 + 半透明误差阴影
%
% 用法：
%   plot_with_shaded_error(x, y_mean, y_err, [R G B])

x = x(:)';
y_mean = y_mean(:)';
y_err = y_err(:)';

y_upper = y_mean + y_err;
y_lower = y_mean - y_err;

x_fill = [x, fliplr(x)];
y_fill = [y_upper, fliplr(y_lower)];

fill(x_fill, y_fill, color, ...
     'FaceAlpha', 0.25, ...
     'EdgeColor', 'none', ...
     'HandleVisibility', 'on');

hold on;
plot(x, y_mean, 'Color', color, 'LineWidth', 5);
end
