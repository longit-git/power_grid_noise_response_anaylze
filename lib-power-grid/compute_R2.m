function R2 = compute_R2(y, yhat)
    % y    : true data (n x 1 or 1 x n)
    % yhat : model prediction

    y = y(:);
    yhat = yhat(:);

    SS_res = sum((y - yhat).^2);         % RSS
    SS_tot = sum((y - mean(y)).^2);      % TSS

    R2 = 1 - SS_res / SS_tot;
end
