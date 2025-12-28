%noise and response for b=1.67 (5/3)
clear;
tic;
graph_name='G20';
N=100000;
h=0.01;
Fs=1/h;
FQ=Fs/N:(Fs/N):Fs;
load(['./currently-using-key-data/',graph_name,'.mat'])
edgelist=G.Edges{:,1:2};
edgelist=edgelist(:,1:2);
load('./currently-using-key-data/pos_noise.mat');
 
num_nodes=max(max(edgelist));
num_edges=length(edgelist);


b=1.67;
w=0.1;
noise=cn_machine_uniform(2*N,b,w,Fs);
noise=(0.3/std(noise))*noise; %0.3为数值模拟中P_i的30%
c=var(noise)/sum((FQ.^(-b))./2);
params_eqs=struct('alpha',1);
params_misc_numerical=struct('time_step',h,'num_nodes',num_nodes,...
    'num_edges',num_edges,'pos_noise',pos_noise,'FQ',FQ);
[t,y_der,evec,eval]=numerical_Response_no_interp(edgelist,noise,pos_noise,params_eqs,params_misc_numerical);
fn=sprintf(['./data/data_of_',graph_name,'/1p67.mat']);
[dirname, ~, ~]=fileparts(fn);
mkdir(dirname);
save(fn,'y_der')
toc;
writematrix(y_der,['./data/data_of_',graph_name,'/1p67.dat'],'Delimiter','\t');
