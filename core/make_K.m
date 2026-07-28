function K = make_K(edgelist, num_nodes, num_edges, mode)
% MAKE_K  生成耦合强度矩阵（num_nodes×num_nodes，对称）
%
% mode:
%   'fixed_30'   — 所有边 K=30
%   'rand_25_35' — 每条边 K~U[25,35]
%   'rand_20_40' — 每条边 K~U[20,40]

K = zeros(num_nodes, num_nodes);

switch mode
    case 'fixed_30'
        for e = 1:num_edges
            i = edgelist(e,1); j = edgelist(e,2);
            K(i,j) = 30;  K(j,i) = 30;
        end

    case 'rand_25_35'
        for e = 1:num_edges
            i = edgelist(e,1); j = edgelist(e,2);
            Kr = 25 + 10*rand();
            K(i,j) = Kr;  K(j,i) = Kr;
        end

    case 'rand_20_40'
        for e = 1:num_edges
            i = edgelist(e,1); j = edgelist(e,2);
            Kr = 20 + 20*rand();
            K(i,j) = Kr;  K(j,i) = Kr;
        end

    otherwise
        error('make_K: 未知 mode="%s"', mode);
end
end
