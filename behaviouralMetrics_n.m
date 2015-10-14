function metric = behaviouralMetrics_n(data, func, nobs)
% function behaviouralMetrics_n(data, func, nobs)
% Runs behavioralMetrics(data, func) with nobs randomly sampled (without
% replacement) observations from data. 

    NOBS = size(data,1);
    if nobs > NOBS; nobs = NOBS; end

    i = randsample(NOBS, nobs);
    metric = behaviouralMetrics(data(i,:), func);

end


