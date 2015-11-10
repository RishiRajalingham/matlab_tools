function [rs_mu, rs_sig, rs_dist] = splitHalfInternalConsistency(data, func, splits_M, splits_type)
% function [rs_mu, rs_sig] = splitHalfInternalConsistency(data, func, M)
% data is standard [s,m,nm,sel], func flag, 
% M: numTrials/condition to consider, max possible if undefined

    corrtype = getCorrelationType();
    if isempty(corrtype); corrtype = 'spearman'; end;
    N = getNumSplitHalves();
    if isempty(N); N = 100; end;

    
    s = data(:,1); % sample category
    m = data(:,2); % match test category
    nm = data(:,3); % nonmatch test category
    us = unique(s); % all categories
    Ns = length(us); % number of categories
    
    if size(data,1) > 10^6; N = 100; end;
    
    if ~exist('splits_M', 'var') || isempty(splits_m)
        splits_M = 10^6; % max value of M;
        splits_type = [];
    end
    % use precomputed splits or recompute them
    if isa(splits_M, 'double')
        M = splits_M;
        splits = getSplitHalves(data, N, M, splits_type);
    elseif isa(splits_M, 'cell')
        splits = splits_M;
        N = length(splits);
    end
    
    rs_dist = zeros(N,1);
    for ni = 1:N
        rs_dist(ni) = getInternalConsistency(splits{ni});
    end

    [rs_mu, rs_sig] = grpstats(rs_dist, [], {'mean', 'std'});
    
    % compute metrics on splits and correlate
    
    function rs = getInternalConsistency(ind)
        
        metric = [];
        for i = 1:2
            ii = ind{i};
            metric = [metric, behaviouralMetrics(data(ii,:), func)];
        end
    
        t = isfinite(metric(:,1)) & isfinite(metric(:,2)); 
        if sum(t) == 0; 
            display('Splithalf error: Not enough data.'); 
            rs = nan;
        else
            rs = corr(metric(t,:), 'type', corrtype); 
            rs = sb_correct(rs(1,2));
        end
    end
    
    % spearman brown correction
    function rho2 = sb_correct(rho)
        rho2 = 2*rho / (1+rho);
    end
    
end