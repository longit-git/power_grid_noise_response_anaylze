% Prepare data by running "do_pd_storage.m"
% assess heterogeneity by Hellinger distance.
graph_name='G20';
load(['./data/data_of_',graph_name,'/pd_storage/pd_storage.mat'])
HD=zeros(21,num_nodes);
for i_b=1:21
    for i_n=1:num_nodes
        ACC1=0;
        for i_k=1:num_realizations-1
            for j=i_k+1:num_realizations
                P1=squeeze(PD_storage(i_b,i_k,i_n,2,:));
                Q1=squeeze(PD_storage(i_b,j,i_n,2,:));
                ACC1=ACC1+Hellinger_distance(P1,Q1);
            end
        end
        HD(i_b,i_n)=ACC1/((num_realizations)*(num_realizations-1)/2);
    end
    disp(['stage3:',num2str(i_b)]);
end
save(['./small_data_for_plotting/',graph_name,'/HD.mat'],'HD');
