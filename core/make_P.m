function P = make_P(num_nodes, mode)
% MAKE_P  生成节点功率分配向量（列向量，长度 num_nodes，总功率守恒=0）
%
% mode:
%   'default' — +1/-1 交替，各占50%；奇数节点时末节点P=0
%   '20_80'   — 20% 节点 P=+0.8，80% 节点 P=-0.2；要求 num_nodes 整除5
%   '40_60'   — 40% 节点 P=+0.6，60% 节点 P=-0.4；要求 num_nodes 整除5

switch mode
    case 'default'
        P = zeros(num_nodes, 1);
        for i = 1 : floor(num_nodes/2)
            P(2*i-1) = +1;
            P(2*i)   = -1;
        end

    case '20_80'
        if mod(num_nodes, 5) ~= 0
            error('make_P: 20_80 模式要求 num_nodes 能被5整除，当前=%d', num_nodes);
        end
        n_pos = round(0.2 * num_nodes);
        P = zeros(num_nodes, 1);
        P(1:n_pos)           = +0.8;
        P(n_pos+1:num_nodes) = -0.2;

    case '40_60'
        if mod(num_nodes, 5) ~= 0
            error('make_P: 40_60 模式要求 num_nodes 能被5整除，当前=%d', num_nodes);
        end
        n_pos = round(0.4 * num_nodes);
        P = zeros(num_nodes, 1);
        P(1:n_pos)           = +0.6;
        P(n_pos+1:num_nodes) = -0.4;

    otherwise
        error('make_P: 未知 mode="%s"', mode);
end
end
