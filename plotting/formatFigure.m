function formatFigure(ax, fontSize, lineWidth, figSize, axPos)
% formatFigure 统一设置图像字体大小、边框和数轴厚度
%
% 用法：
%   formatFigure(gca, 14, 1.5);
%
% 输入参数：
%   ax        - 坐标轴句柄 (一般用 gca)
%   fontSize  - 字体大小 (默认 12)
%   lineWidth - 轴线与边框厚度 (默认 1.2)

    if nargin < 2 || isempty(fontSize)
        fontSize = 20;
    end
    if nargin < 3 || isempty(lineWidth)
        lineWidth = 1.2;
    end

    % 设置字体
    set(ax, 'FontSize', fontSize, ...
            'FontName', 'Helvetica', ... % 可改为 Arial/Helvetica
            'LineWidth', lineWidth, ...
            'Box', 'on'); % 开启边框

    % 坐标轴标签字体大小
    ax.XLabel.FontSize = fontSize;
    ax.YLabel.FontSize = fontSize;
    ax.Title.FontSize  = fontSize + 2;

    % 如果有图例，也设置
    leg = findobj(ax.Parent, 'Type', 'Legend');
    if ~isempty(leg)
        set(leg, 'FontSize', fontSize, 'LineWidth', lineWidth);
    end

    % 如果有颜色条，也设置
    cb = findobj(ax.Parent, 'Type', 'ColorBar');
    if ~isempty(cb)
        set(cb, 'FontSize', fontSize, 'LineWidth', lineWidth);
    end

    % 控制图像大小
    if nargin < 4 || isempty(figSize)
        figSize=[400, 300];
    end
        fig = ancestor(ax, 'figure');
        oldPos = fig.Position;  % [x0 y0 width height]
        fig.Position = [oldPos(1), oldPos(2), figSize(1), figSize(2)];
    % 控制图像大小
    if nargin < 5 || isempty(axPos)
        axPos=[0.15, 0.25, 0.68, 0.68];
    end
    ax.Position = axPos;
end

