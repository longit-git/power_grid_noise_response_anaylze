function output_parameter(fid,graph_name,N,h,b,num_realizations,pos_noise,interp_or_not)
fprintf(fid,['graph_name=',graph_name,'\n']);
fprintf(fid,'N=%d\n',N);
fprintf(fid,'h=%.2f\n',h);
fprintf(fid,'b=%.2f\n',b);
fprintf(fid,'times of realization=%d\n',num_realizations);
fprintf(fid,'pos_noise=%d\n',pos_noise);
fprintf(fid,'#if pos_noise=1, as it should be, the perturbed node is the "first node" in network generating process\n');
fprintf(fid,interp_or_not);
end
