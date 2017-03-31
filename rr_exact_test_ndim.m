function [p,z,dist_data] = rr_exact_test_ndim(x_null, x_obs, dist_met)
% function [p,z,dist_data] = rr_exact_test_ndim(x_null, x_obs, dist_met)
% multidimensional exact test: for a given multidimensional distribution 
% x_null, compute empirical probability of observing x_obs, based on 
% proportion of **heldout** samples that are at a greater distance.
% Distance can be computed using Mahalanobis and euclidean metrics.
% Null distributions are fitted to functional distributions (normal, gamma)
% and the best fitting function (based on KL divergence) is used to compute
% a p value.

    [ns, nd] = size(x_null);
    if ~exist('x_obs', 'var') || isempty(x_obs)
        x_obs = zeros(1,nd);
    end
    if ~exist('dist_met', 'var')
        dist_met = 'mahal';
    end
    
    assert(size(x_null,2) == length(x_obs));
    if size(x_obs,1) ~= 1
        x_obs = x_obs';
    end

    % remove nans
    t = ~isnan(x_obs);
    x_null = x_null(:,t);
    x_obs = x_obs(:,t);
    
    % estimate mahal distance using parameters (mu,sig) from different
    % samples.
    x_train = x_null(1:floor(ns/2),:);
    
    x_test = x_null(1+floor(ns/2):end,:);
    ntrain = size(x_train,1);
    ntest = size(x_test,1);
    d_dist = nan(ntest,1);
    
    % remove dimensions with no variance
    dim_oi = var(x_train) > 0;
    x_train = x_train(:, dim_oi);
    x_test = x_test(:, dim_oi);
    x_obs = x_obs(:, dim_oi);
    
    if length(x_obs) == 1
        d_0 = x_obs;
        d_dist = cat(1, x_train, x_test);
        pd = determine_distribution(d_dist, {'Normal'});
    else
        if strcmp(dist_met, 'mahal')
            d_0 = sqrt(mahal(x_obs, x_train)); %mahalanobis distance^2 from origin
            for i = 1:ntest
                d_dist(i) = sqrt(mahal(x_test(i,:), x_train));
            end
            pd = determine_distribution(d_dist, {'Gamma', 'Normal'});
        elseif strcmp(dist_met, 'euclid')
            D = squareform(pdist([x_obs; x_test; x_train]));
            d_0 = nanmean(D(1,ntest+2:end),2);
            d_dist = nanmean(D(2:(ntest+1),ntest+2:end),2);
            pd = determine_distribution(d_dist, {'Gamma', 'Normal'});
        end
    end
    
    del_dist = repmat(d_0,size(d_dist,1),1) - d_dist;
    dist_data.d_dist = d_dist;
    dist_data.d_0 = d_0;
    dist_data.del_dist = del_dist;
    p = 1 - cdf(pd, d_0);
    z = (d_0 - nanmean(d_dist)) ./ nanstd(d_dist);
        
    function pd = determine_distribution(d, dist_opts)
        if ~exist('dist_opts', 'var')
            dist_opts = {'Normal', 'Gamma'};
        end
        xvals = linspace(min(d), max(d), 10);
        emp_dist = histc(d, xvals);
        emp_dist = emp_dist ./ sum(emp_dist);
        
        
        ND = length(dist_opts);
        PDs = cell(ND,1);
        KL_score = nan(ND,1);
        for doi = 1:length(dist_opts)
            PDs{doi} = fitdist(d, dist_opts{doi});
            fit_dist = pdf(PDs{doi}, xvals);
            KL_score(doi) = KLDiv(fit_dist(:)', emp_dist(:)');
        end
        
        [~,mi] = min(KL_score);
        pd = PDs{mi};
    end
    
end