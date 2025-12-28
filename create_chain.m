%% chain network
N = 10;

% 构造链式邻接矩阵
A = diag(ones(N-1,1),1) + diag(ones(N-1,1),-1);

% 创建图
G = graph(A);

% 绘制
figure;
plot(G,'Layout','force');
title('链式网络');

