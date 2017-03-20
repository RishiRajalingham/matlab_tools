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

    % estimate mahal distance using parameters (mu,sig) from different
    % samples.
    x_train = x_null(1:floor(ns/2),:);
    x_test = x_null(1+floor(ns/2):end,:);
    ntrain = size(x_train,1);
    ntest = size(x_test,1);
    d_dist = nan(ntest,1);

    if strcmp(dist_met, 'mahal')
        d_0 = sqrt(mahal(x_obs, x_train)); %mahalanobis distance^2 from origin
        for i = 1:ntest
            d_dist(i) = sqrt(mahal(x_test(i,:), x_train));
        end
    elseif strcmp(dist_met, 'euclid')
        D = squareform(pdist([x_obs; x_test; x_train]));
        d_0 = nanmean(D(1,end-ntrain+1:end));
        d_dist = nanmean(D(2:ntest+1,end-ntrain+1:end))';
    end

    del_dist = repmat(d_0,ntest,1) - d_dist;
    dist_data.d_dist = d_dist;
    dist_data.d_0 = d_0;
    dist_data.del_dist = del_dist;
    
    pd = determine_distribution(d_dist);
    p = 1 - cdf(pd, d_0);
    z = (d_0 - nanmean(d_dist)) ./ nanstd(d_dist);
    
    
%     check_normality = ~kstest(zscore(del_dist));
%     if check_normality
%         dist_to_use = 'Normal';
%     end

    
    
    
    function pd = determine_distribution(d)
        xvals = linspace(min(d), max(d), 50);
        emp_dist = histc(d, xvals);
        emp_dist = emp_dist ./ sum(emp_dist);
        
        dist_opts = {'Normal', 'Gamma'};
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