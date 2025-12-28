% prepare data by running do_pd_storage.m
clear
graph_name='G20';
load(['./data/data_of_',graph_name,'/pd_storage/pd_storage.mat'])
HD_vs_num_realizations=zeros(21,num_nodes,num_realizations);
for r_length=2:1:num_realizations
    for i_b=1:21
        for i_n=1:num_nodes
            ACC1=0;
            for i_k=1:r_length-1
                for j=i_k+1:r_length
                    P1=squeeze(PD_storage(i_b,i_k,i_n,2,:));
                    Q1=squeeze(PD_storage(i_b,j,i_n,2,:));
                    ACC1=ACC1+Hellinger_distance(P1,Q1);
                end
            end
            HD_vs_num_realizations(i_b,i_n,r_length)=ACC1/((r_length)*(r_length-1)/2);
        end
    end
    disp(['stage3:',num2str(r_length)]);
end

%%
fn_HD=['small_data_for_plotting/',graph_name,'/HD_vs_num_realizations.mat'];
[fn_dir,~,~]=fileparts(fn_HD);
mkdir(fn_dir);
save(['small_data_for_plotting/',graph_name,'/HD_vs_num_realizations.mat'],'HD_vs_num_realizations');


