function result=Hellinger_distance(P,Q)
% Input P and Q are two sets of possibility distribution
result=sqrt(sum((sqrt(P)-sqrt(Q)).^2))/(sqrt(2));
end