% This script stores possibility density of responses on every node, with
% every b, every realization
% to be used to calculate "H".
clear
num_bin = 35;
num_realizations=30;
num_nodes = 200;
num_b = 21;
graph_name='G20';
PD_storage = zeros(num_b, num_realizations, num_nodes, 2, num_bin);
edge_limit = zeros(num_b, num_nodes, 2);

max1 = zeros(num_realizations, num_nodes);
min1 = zeros(num_realizations, num_nodes);

%% ---------------------------
%   Stage 1: Compute edge limits
% ----------------------------
for i_b = 1:num_b
    for i_k = 1:num_realizations

        fn= sprintf(['./data/data_of_',graph_name,'/y_der_%d_%d.mat'],i_b,i_k);
        data = load(fn);
        X = data.y_der_w;

        max1(i_k, :) = max(X, [], 1);
        min1(i_k, :) = min(X, [], 1);
    end

    edge_limit(i_b, :, 1) = min(min1, [], 1);
    edge_limit(i_b, :, 2) = max(max1, [], 1);

    fprintf('Stage1 finished: %d\n', i_b);
end

%% ---------------------------
%   Stage 2: Compute histograms
% ----------------------------
parfor i_b = 1:num_b
    local_PD = zeros(num_realizations, num_nodes, 2, num_bin);

    for i_k = 1:num_realizations
        fn= sprintf(['./data/data_of_',graph_name,'/y_der_%d_%d.mat'],i_b,i_k);
        X = load(fn).y_der_w;

        for i_n = 1:num_nodes
            edges = linspace(edge_limit(i_b,i_n,1), edge_limit(i_b,i_n,2), num_bin+1);
            nums = histcounts(X(:,i_n), edges, 'Normalization','probability');

            local_PD(i_k, i_n, 1, :) = (edges(1:end-1)+edges(2:end))/2;
            local_PD(i_k, i_n, 2, :) = nums;
        end
    end

    PD_storage(i_b, :, :, :, :) = local_PD;

    fprintf('Stage2 finished: %d\n', i_b);
end

fn = ['./data/data_of_',graph_name,'/pd_storage/pd_storage.mat'];
[fndir,~,~]=fileparts(fn);
mkdir(fndir);
save(fn);