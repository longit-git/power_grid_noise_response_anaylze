function [chi2_stat, r2, rmse] = Gaussianity_chi2(data, num_bin)
% chi2 = sum( (model - y).^2 ./ model )
% model: Gaussian pdf
% y: empirical pdf (histogram)

    [~, column] = size(data);

    chi2_stat = nan(1, column);
    sigma_hat = nan(1, column);
    r2        = nan(1, column);
    rmse      = nan(1, column);

    for fit_node = 1:column
        x = data(:, fit_node);

        % ===== 1. Fit Gaussian =====
        pd = fitdist(x, 'Normal');
        mu    = pd.mu;
        sigma = pd.sigma;
        sigma_hat(fit_node) = sigma;

        % ===== 2. Empirical pdf =====
        [y, edges] = histcounts(x, num_bin, 'Normalization', 'pdf');
        bin_centers = (edges(1:end-1) + edges(2:end)) / 2;

        % ===== 3. Gaussian model (pdf) =====
        model = pdf(pd, bin_centers);

        % ===== 4. valid bins =====
        valid = model > 0;
        y     = y(valid);
        model = model(valid);

        n = numel(y);

        % ===== 5. χ² / R² / RMSE =====
        y_mean = mean(y);

        ss_res = 0.0;
        ss_tot = 0.0;
        chi2   = 0.0;

        for i = 1:n
            residual = model(i) - y(i);
            ss_res = ss_res + residual^2;
            ss_tot = ss_tot + (y(i) - y_mean)^2;
            chi2   = chi2   + residual^2 / model(i);
        end

        r2(fit_node)   = 1.0 - ss_res / ss_tot;
        rmse(fit_node) = sqrt(ss_res / n);
        chi2_stat(fit_node) = chi2;
    end
end
