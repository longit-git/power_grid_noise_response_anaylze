% no interpolation, use real noise data (in line with the "color") to fit
% rk4
function [t,response_y_der,evec,eval]=numerical_Response_no_interp(edgelist,colored_b,pos_noise,paramsEqs,paramsMiscNumerical)

% initialize some parameters
alpha=paramsEqs.alpha;
h=paramsMiscNumerical.time_step;
num_nodes=paramsMiscNumerical.num_nodes;
num_edges=paramsMiscNumerical.num_edges;
P=PDistribution_Costomized(num_nodes);

K=KDistribution(edgelist,num_nodes,num_edges);
N = length(colored_b)/2; % length of times series

y0=zeros(1,2*num_nodes);
system_func=@(array) func_generator(K,P,alpha,array);
y0=fsolve(system_func,y0);

tspan=[0,N*h-h];

% the noise applied on P_i
P_p=@(t) apply_noise(pos_noise,colored_b,P,h,t);
% handle of the model
power_model_p=@(t,y) func_generator(K,P_p(t),alpha,y);
% rk4
[t,y]=rungeKuttaSolver(power_model_p,tspan,y0,h);

[evec,eval]=eig(Lap_Matrix(y0(num_nodes+1:end),K));

response_y_der=y(1:num_nodes,:).';
end

%% function part
% this func represents the powergrid system
% ODE_array, the first half of it represents the derivatives of theta_i, the
% second half represents theta_i.
function ODE_array=func_generator(K,P,alpha,y)
ODE_array=zeros(length(y),1);
for i=1:length(P)
    ODE_array(i)=P(i)-alpha*y(i)+sum_function(i);
end
for i=1:length(P)
    ODE_array(length(P)+i)=y(i);
end
    function s=sum_function(i)
        s=0;
        for j=1:length(P)
            s=s+K(i,j)*sin(y(length(P)+j)-y(length(P)+i));
        end
    end
end

%white Gaussian noise
function array=apply_noise(position_p,noise,P,h,t)
index=round(t/(h/2)+1);
P(position_p)=P(position_p)+noise(index);
array=P;
end
%laplacian matrix
function LapMatrix = Lap_Matrix(theta, K)
num_points = length(theta);
LapMatrix = zeros(num_points, num_points);
for i = 1:num_points
    for j = 1:num_points
        if i == j
            LapMatrix(i, j) = 0;
            for n = 1:num_points
                LapMatrix(i, j) = LapMatrix(i, j) + K(i, n) * cos(theta(n) - theta(i));
            end
        else
            LapMatrix(i, j) = -K(i, j) * cos(theta(j) - theta(i));
        end
    end
end
end
%allocate P
function P_Costomized=PDistribution_Costomized(num_nodes)
if mod(num_nodes,2)==0
    P_Costomized=zeros(num_nodes,1);
    for i=1:num_nodes/2
        P_Costomized(2*i)=-1;
        P_Costomized(2*i-1)=1;
    end
else
    P_Costomized=zeros(num_nodes,1);
    for i=1:floor(num_nodes/2)
        P_Costomized(2*i)=-1;
        P_Costomized(2*i-1)=1;
    end
    P_Costomized(end)=0;
end
end
%K_ij
function K=KDistribution(edgelist,num_nodes,num_edges)
K=zeros(num_nodes,num_nodes);
for unique_i=1:num_edges
    K(edgelist(unique_i,1),edgelist(unique_i,2))=30;
    K(edgelist(unique_i,2),edgelist(unique_i,1))=30;
end
end