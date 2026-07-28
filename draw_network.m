% fisrt, load a Graph object you want to draw; second, run the script.
figure;
p=plot(G,'Layout','force');
regular_node_color=zeros(numnodes(G),3);
regular_node_color(:,1)=0.5;
regular_node_color(:,2)=0.5;
regular_node_color(:,3)=0.5; %color默认值
%regular_node_color(36,:)=[0.48 0.53 0.78];
p.NodeColor=regular_node_color;
load('./currently-using-key-data/pos_noise.mat');
markersize=8*ones(numnodes(G),1);
markersize([pos_noise])=10;
p.MarkerSize=markersize;
p.NodeLabel={};
p.EdgeColor=[0.5,0.5,0.5];
