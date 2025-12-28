function result=Hellinger_distance(P,Q)
result=(sum((sqrt(P)-sqrt(Q)).^2))/(sqrt(2));
end