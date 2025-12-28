% customized graph, no interpolation
clear;
prompt = "!!!Are you sure to cover the previous data?!!! Y/N ";
txt = input(prompt,"s");
if txt == 'Y'
else
    error('Data writting was not permitted');
end

tic;

N=100000;
h=0.01;
Fs=1/h;
FQ=Fs/N:(Fs/N):Fs;
load('./currently-using-key-data/pos_noise.mat');

b=0:0.1:2;
num_realizations=30;
graph_list = {'G20','G21','G_BA','G_chain_10_nodes'};
for i_G=1:length(graph_list)
    graph_name=graph_list{i_G};
    fn_G=sprintf(['./currently-using-key-data/',graph_name,'.mat']);
    load(fn_G);
    edgelist=G.Edges{:,1:2};
    edgelist=edgelist(:,1:2);
    num_nodes=max(max(edgelist));
    num_edges=length(edgelist);
    mkdir(['./data/data_of_',graph_name]);
    for i_b=1:length(b)
        noise_rtw=zeros(2*N,num_realizations);
        y_der_rtw=zeros(N,num_nodes,num_realizations);
        parfor i_k=1:num_realizations
            w=0.1;
            noise=cn_machine_uniform(2*N,b(i_b),w,Fs);
            noise=(0.3/std(noise))*noise; %0.3为数值模拟中P_i的30%
            noise_rtw(:,i_k)=noise;
            c=var(noise)/sum((FQ.^(-b(i_b)))./2);
            params_eqs=struct('alpha',1);
            params_misc_numerical=struct('time_step',h,'num_nodes',num_nodes,...
                'num_edges',num_edges,'pos_noise',pos_noise,'FQ',FQ);
            [t,y_der,evec,eval]=numerical_Response_no_interp(edgelist,noise,pos_noise,params_eqs,params_misc_numerical);
            y_der_rtw(:,:,i_k)=y_der;
            params_misc_numerical.evec=evec;
            params_misc_numerical.eval=eval;
        end
        for i_write=1:num_realizations
            y_der_w=y_der_rtw(:,:,i_write);
            noise_w=noise_rtw(:,i_write);
            fname_y_der=sprintf(['./data/data_of_',graph_name,'/y_der_%d_%d.mat'],i_b,i_write);
            save(fname_y_der,'y_der_w');
            fname_noise=sprintf(['./data/data_of_',graph_name,'/noise_%d_%d.mat'],i_b,i_write);
            save(fname_noise,'noise_w')
        end
        disp(graph_name);
        disp(num2str(i_b));
        fid=fopen(['./data/data_of_',graph_name,'/parameter.txt'],'w');
    end
    output_parameter(fid,graph_name,N,h,b,num_realizations,pos_noise,'no interpolation');
end
toc;
