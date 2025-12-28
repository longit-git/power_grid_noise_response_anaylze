
function [non_Gaussianity,width]=Gaussianity_assessment(data,num_bin)
    % input
    % data: n * column
    % num_bin: number of bins
    sizeData=size(data);
    row=sizeData(1);
    column=sizeData(2);
    non_Gaussianity=zeros(1,column);
    width=zeros(1,column);
    for fit_node=1:column
    sub_data = data(:,fit_node);
    % fit Gaussian distribution
    pd = fitdist(sub_data, 'Normal');
    % PDF
    [bin_counts, bin_edges] = histcounts(sub_data, num_bin, 'Normalization', 'pdf');
    
    % center of bins
    
    bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;
    bin_width=bin_centers(2)-bin_centers(1);
    hold on;

    % RMSE
    err=bin_width*rmse(bin_counts(2:end-1),pdf(pd,bin_centers(2:end-1)));
    non_Gaussianity(:,fit_node)=err^2;
    width(:,fit_node)=pd.sigma;
    end
end