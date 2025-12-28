% create an random power grid
% method from: Schultz, P., Heitzig, J. & Kurths, J. A random growth model 
% for power grids and other spatially embedded infrastructure networks. 
% Eur. Phys. J. Spec. Top. 223, 2593–2610 (2014). 
% https://doi.org/10.1140/epjst/e2014-02279-6
clear
N = 200;
N0 = 5;
p = 0.8;
q = 0.8;
r = 28;
s = 0.8;
G = init(N0,p,q,r,s);
G = grow(G,N-N0,p,q,r,s);

%print statistics(actual mean degree, actual path length, transitivity(/mean degree),ratio r{k>kbar},degree correlation,Resulting gamma from linear fit to log-decumulative,Fiedler eigenvalue lambda2)

%计算平均度
mean_degree = mean(degree(G));
plot(G)

 
%以下为函数
function point = uniform_unitsquare()
%return point drawn uniformly at random from unit square
point = rand(2,1);
end
%
function dist = euclidean(x,y)
% return Euclidean distance between x and y
dist = sqrt(sum((x-y).^2));
end

% model components:
function G_out = Graph_full(N0)
% Return complete graph on N0 nodes
G = graph();
for i = 1:N0
	for j = i+1:N0
		G = addedge(G,i,j);
	end
end
G_out = G;
end
%接下来定义相应的初始化模块,生长模块
function G_out =  init(N0,p,q,r,s)
rho = @uniform_unitsquare;
d = @euclidean;
%%
%    return initialized network (an igraph Graph with additional attributes like coordinates)

%    :param N0: initial number of nodes
%    :param p: attachment probability [0;1]
%    :param q: attachment probability [0;1]
%    :param r: distance redundancy trade-off exponent
%    :param s: edge splitting probability [0;1]
%    :param rho: distribution to draw node locations from (function argument)
%    :param d: distance metric (function argument)
%    :return: G: Graph object
%%

%第一步:使用rho生成N0x2的矩阵,每列代表一个节点的位置
x = zeros(N0,2);
for i = 1:N0
	x(i,:) = rho();%结果:x是N0x2的矩阵
end
%第二步:构建最小生成树
full_graph = Graph_full(N0);
factor = 1e5;

%计算每条边的权重
weights = zeros(numedges(full_graph),1);
for edge_count = 1:numedges(full_graph)
	edge = full_graph.Edges(edge_count, :);
	source = edge.EndNodes(1);
	target = edge.EndNodes(2);
	weights(edge_count) = factor * d(x(source, :), x(target, :));
end
%对图G赋予权重
full_graph.Edges.Weights =weights;
G = minspantree(full_graph);
%删除中间变量
clear full_graph;
%第三步: add redundant links
m = min(floor(N0*(1 - s)*(p + q)), floor(N0*(N0 - 1)/2 - (N0 - 1)));
for i = 1:N0
	for j = 1:N0
		G.Nodes.dmatrix(i,j) = d(x(i,:),x(j,:));
	end
end
for i = 1:N0
	G.Nodes.dmatrix(i,i) = inf;
end
for edge_count = 1:numedges(G)
	edge = G.Edges(edge_count,:);
	i = edge.EndNodes(1);
	j = edge.EndNodes(2);
	G.Nodes.dmatrix(i,j) = inf;
	G.Nodes.dmatrix(j,i) = inf;
end
if r>0
	G.Nodes.dGmatrix = distances(G,'Method','unweighted');
end

for k = 1:m
	if r > 0
		target_matrix = G.Nodes.dmatrix./(G.Nodes.dGmatrix+1).^r;
	else
		target_matrix = G.Nodes.dmatrix;
	end

	%找到target_matrix上最小值对应的i,j索引,添加到G
	[i,j] = find(target_matrix == min(target_matrix(:)));
	i = i(1);
	j = j(1);
	G = addedge(G,i,j);
	G.Nodes.dmatrix(i,j) = inf;
	G.Nodes.dmatrix(j,i) = inf;
	%调整网络距离
	if r > 0
		G.Nodes.dGmatrix = min(min(G.Nodes.dGmatrix,G.Nodes.dGmatrix(:, i)+ones(N0, N0)+G.Nodes.dGmatrix(j, :)),G.Nodes.dGmatrix(:, j) + ones(N0, N0) + G.Nodes.dGmatrix(i, :));
	end

	if r > 0
		assert(all(all(abs(G.Nodes.dGmatrix - distances(G,'Method','unweighted')) < 1e-10)), 'Check that update went well');
	end
end
%返回值
G.Nodes.x = x;
clear x;
G_out = G;
end

%%

function G_out = grow(G,n,p,q,r,s)
% return graph G with n additional nodes(initially, there are N, so you end up with N+n nodes)
%G:当前网络(Graph object)
%n:增加的节点数
%p:
%q:
%r:distance / redundancy trade-off exponent
%s: 边分裂概率
%rho: 位置生成函数
%d: 距离度规
%返回值:增加节点后的网络(Graph object)

rho = @uniform_unitsquare;
d = @euclidean;
x = zeros(n,2);
for i = 1:n
	x(i,:) = rho();
end
%更新图G的节点数
%准备新的的网络数据
%获取图G的节点数
N = numnodes(G);
G = addnode(G,n);
G.Nodes.dmatrix = [G.Nodes.dmatrix(1:N,1:N),inf(N,n);inf(n,N+n)];%这里一定事先要把dmatrix的大小限制到NxN,因为添加节点后,MATLAB会自动扩展矩阵行数,但不会自动扩展列数
G.Nodes.x((N+1):(N+n),:) = x;
if r>0
	G.Nodes.dGmatrix = [G.Nodes.dGmatrix(1:N,1:N),zeros(N,n);zeros(n,N+n)];
end
onesquare = ones(N+n,N+n);
counter = 1;
for i = (N+1):(N+n)
	counter = counter + 1;
	if (rand(1) < s) && (isempty(G.Edges) == 0)
		%step G5 split random edge at midpoint

		%choose random edge
		ei = randi(height(G.Edges));
		e = G.Edges(ei,:);
		%get the two nodes
		a = e.EndNodes(1);
		b = e.EndNodes(2);

		%add node at midpoint and calc distances:

		G.Nodes.x(i,:) = (G.Nodes.x(a,:)+G.Nodes.x(b,:))/2;
		G.Nodes.dmatrix(i,1:N) = 1; %需要考证
		G.Nodes.dmatrix(1:N,i) = 1; %需要考证

		%replace edge with two edges
		G = rmedge(G,a,b);% Remark: one of the lists G.xx is not updated
		G = addedge(G,a,i);
		G = addedge(G,i,b);

		%make sure i,a and i,b are not selected again
		G.Nodes.dmatrix(i,a) = inf;
		G.Nodes.dmatrix(a,i) = inf;
		G.Nodes.dmatrix(i,b) = inf;
		G.Nodes.dmatrix(b,i) = inf;

		%recalculate all network distances:
		if r>0
			G.Nodes.dGmatrix = distances(G,'Method','unweighted');%需要考证
		end
	else
		%step G2: link to nearest node

		% 计算第 i 行的前 N 列
		G.Nodes.dmatrix(i, 1:N) = arrayfun(@(j) d(G.Nodes.x(i), G.Nodes.x(j)), 1:N);

		% 计算前 N 行的第 i 列
		G.Nodes.dmatrix(1:N, i) = arrayfun(@(j) d(G.Nodes.x(i), G.Nodes.x(j)), 1:N);%需要考证

		di = G.Nodes.dmatrix(i,1:N);
		[~, j] = min(di);
		G = addedge(G,i,j);
		%make sure i,j are not selected again:
		di(j) = inf;
		G.Nodes.dmatrix(i,j) = inf;
		G.Nodes.dmatrix(j,i) = inf;
		%adjust network distances:
		if r>0
			G.Nodes.dGmatrix(i,1:N) = G.Nodes.dGmatrix(j,1:N) + 1;
			G.Nodes.dGmatrix(1:N,i) = G.Nodes.dGmatrix(j,1:N) + 1;
			dGi = G.Nodes.dGmatrix(i,1:N);
		end


		%step G3: add optimal redundant edge to i
		if rand(1) < p
			if r >0
				target_vector = di./(dGi+1).^r;
			else
				target_vector = di;
			end
			%find minimum and add to G:
			[~,k] = min(target_vector);
			if target_vector(k) < inf % otherwise i is already linked to all other nodes
				assert(~findedge(G,i,k));%待定
				G = addedge(G,i,k);
				%make sure i,k are not selected again:
				G.Nodes.dmatrix(k,i) = inf;
				G.Nodes.dmatrix(i,k) = inf;
				%adjust network distances:
				if r > 0
					G.Nodes.dGmatrix = min(min(G.Nodes.dGmatrix,G.Nodes.dGmatrix(:, i) + onesquare + G.Nodes.dGmatrix(k, :)),G.Nodes.dGmatrix(:, k) + onesquare + G.Nodes.dGmatrix(i, :));
				end
			end
			counter = counter + 1;
		end


		%step G4 add another optimal redundant link to random node
		if rand(1)<q
			i2 = randi(N);
			if r>0
				target_vector = G.Nodes.dmatrix(i2,1:N)./(G.Nodes.dGmatrix(i2,1:N)+1).^r;
			else
				target_vector = G.Nodes.dmatrix(i2,1:N);
			end
			%find minimum and add to G:
			[~,k2] = min(target_vector);
			if target_vector(k2)<inf%otherwise i2 is already linked to all other nodes
				assert(~findedge(G,i2,k2))
				G = addedge(G,i2,k2);
				G.Nodes.dmatrix(i2,k2)=inf;
				G.Nodes.dmatrix(k2,i2)=inf;
				if r > 0
					G.Nodes.dGmatrix = min(min(G.Nodes.dGmatrix,G.Nodes.dGmatrix(:,i2)+onesquare+G.Nodes.dGmatrix(k2,:)),G.Nodes.dGmatrix(:,k2)+onesquare+G.Nodes.dGmatrix(i2,:));
				end
			end
			counter = counter+1;
		end


	end
	N = N+1;

end
if r > 0
	%assert(all(all(abs(G.Nodes.dGmatrix - distances(G,'Method','unweighted')) < 1e-10)), 'Check that update went well');
end
G_out = G;
end


