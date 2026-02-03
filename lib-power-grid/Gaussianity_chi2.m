function [chi2_stat, sigma_hat, df] = Gaussianity_chi2(data, num_bin)
% data: n x column
% num_bin: number of bins (k)
% chi2_stat: 1 x column chi-square statistic
% sigma_hat: 1 x column fitted sigma
% df: 1 x column degrees of freedom (approx k-3, after filtering bins)

    [n, column] = size(data);
    chi2_stat = nan(1, column);
    sigma_hat = nan(1, column);
    df = nan(1, column);

    for fit_node = 1:column
        x = data(:, fit_node);

        % Fit Normal(mu, sigma)
        pd = fitdist(x, 'Normal');
        sigma_hat(fit_node) = pd.sigma;

        % Histogram counts (NOT pdf)
        [O, edges] = histcounts(x, num_bin, 'Normalization', 'count');

        % Expected counts from fitted CDF over bin edges
        p = cdf(pd, edges(2:end)) - cdf(pd, edges(1:end-1));   % 1 x k
        E = n * p;                                             % 1 x k

        % Remove bins with E=0 (or too small)
        valid = E > 0;   % you can tighten to E >= 5 if you want
        O = O(valid);
        E = E(valid);

        % Pearson chi-square
        chi2 = sum((O - E).^2 ./ E);

        % degrees of freedom: (#valid bins) - 1 - (#estimated params=2)
        k_eff = numel(E);
        df_eff = k_eff - 1 - 2;

        chi2_stat(fit_node) = chi2;
        df(fit_node) = df_eff;
    end
end