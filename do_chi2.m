clear
load('./default_init.mat');
num_b=length_b;
graph_list={'G20','G_chain_10_nodes','G21','G_BA'};
for i_G=1:length(graph_list)
    graph_name=graph_list{i_G};
    load(['./currently-using-key-data/',graph_name,'.mat']);
    num_nodes=numnodes(G);
    chi2=zeros(num_realizations,num_b,num_nodes);
    for i_b=1:num_b
        parfor i_k=1:num_realizations
            fn=sprintf(['./data/data_of_',graph_name,'/y_der_%d_%d.mat'],i_b,i_k);
            y=load(fn);
            X=y.y_der_w;
            chi2(i_k,i_b,:)=Gaussianity_assessment(X,num_bin);
        end
        disp(graph_name);
        disp(i_b);
    end
    fn1=['./data/data_of_',graph_name,'/chi2/chi2.mat'];
    [fn1dir,~,~]=fileparts(fn1);
    mkdir(fn1dir);
    save(fn1,'chi2');

    fn2=['./small_data_for_plotting/',graph_name,'/chi2.mat'];
    [fn2dir,~,~]=fileparts(fn2);
    mkdir(fn2dir);
    save(fn2,'chi2');
end