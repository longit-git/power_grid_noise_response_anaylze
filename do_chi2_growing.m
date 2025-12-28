clear;
graph_name='G20';
load('./default_init.mat');
special_b=[1,6,11,16,21];
num_slides=10;
num_realizations=30;
slide_width=N/10;
slide_ends=slide_width*(1:num_slides);
chi2_growing=zeros(num_slides,num_realizations,length(special_b),num_nodes);
small_figure_node=37;
num_bin=35;
colortable=[0.90, 0.80, 0.40;0.40, 0.70, 0.40;0.19,0.49,0.57];
colortable=[colortable(1,:);(colortable(1,:)+colortable(2,:))/2;colortable(2,:);(colortable(2,:)+colortable(3,:))/2;colortable(3,:)];

%%
for i_b=special_b
    for i_k=1:num_realizations
        fname_noise=sprintf(['./data/data_of_',graph_name,'/noise_%d_%d.mat'],i_b,i_k);
        load(fname_noise);
        fname_y_der=sprintf(['./data/data_of_',graph_name,'/y_der_%d_%d.mat'],i_b,i_k);
        load(fname_y_der);
        if ismember(i_b,special_b)
            for i_slides=1:num_slides
                chi2_growing(i_slides,i_k,(i_b+4)/5,:)=Gaussianity_assessment(y_der_w(1:slide_ends(i_slides),:),num_bin);
            end
        end
        disp('c');
    end
end
%%
fn_growing=['./data/data_of_',graph_name,'/chi2_growing/chi2_growing.mat'];
[fndir,~,~]=fileparts(fn_growing);
mkdir(fndir);
save(fn_growing);



