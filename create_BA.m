%% BA scale-free
% 节点数
N = 200;    
% 初始小网络的节点数
m0 = 4;    
% 每次新增节点的连接数
m = 2;     

% 初始化一个完全图 (m0 节点)
A = ones(m0) - eye(m0);

% 循环添加节点
for newNode = (m0+1):N
    degree = sum(A,2);             % 当前度
    P = degree / sum(degree);      % 连接概率
    
    % 选择 m 个旧节点连接 (按概率)
    targets = randsample(1:(newNode-1), m, true, P);
    
    % 更新邻接矩阵
    A(newNode, targets) = 1;
    A(targets, newNode) = 1;
end

% 创建图
G = graph(A);

% 绘制
figure;
plot(G,'Layout','force');
title('BA无标度网络');
