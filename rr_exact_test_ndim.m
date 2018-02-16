function [p,z,dist_data] = rr_exact_test_ndim(x_null, x_obs, opt)
% function [p,z,dist_data] = rr_exact_test_ndim(x_null, x_obs, dist_met)
% multidimensional exact test: for a given multidimensional distribution 
% x_null, compute empirical probability of observing x_obs, based on 
% proportion of **heldout** samples that are at a greater distance.
% Distance can be computed using any pairwise distance metrics (e.g.
% euclidean, standardized euclidean, mahalanobis, etc...).
% Null distributions are fitted to functional distributions (normal, gamma)
% and the best fitting function (based on KL divergence) is used to compute
% a p value.

    [ns, nd] = size(x_null);
    if ~exist('x_obs', 'var') || isempty(x_obs)
        x_obs = zeros(1,nd);
    end
    
    if ~exist('opt', 'var')
        dist_met = 'euclidean';
        run_pca = 1;
        dist_opts = {'Gamma'};%'{'Normal', 'Gamma'};%, 'tlocationscale'};
    else
        dist_met = opt.dist_met;
        run_pca =  opt.run_pca;
        dist_opts = opt.dist_opts;
    end
    
    assert(size(x_null,2) == length(x_obs));
    if size(x_obs,1) ~= 1
        x_obs = x_obs';
    end

    % remove nans
    t = ~isnan(x_obs);
    x_null = x_null(:,t);
    x_obs = x_obs(:,t);
    
    % estimate  distance from different heldout samples (k-fold).
    % do this niter times on different splits of heldout samples.
    
    niter = 10;
    k = 2;
    nt = ns;
    d_0 = nan(niter,k,nt);
    d_dist = nan(niter,ns,nt);
    for iter = 1:niter
        % randomly split into train test splits
        inds = crossvalind('Kfold',ns,k);
        for ki = 1:k
        
            tr_inds = inds ~= ki;
            te_inds = inds == ki;
            x_train = x_null(tr_inds,:);
            x_test = x_null(te_inds,:);
            x_ob = x_obs;
            ntrain = size(x_train,1);
            ntest = size(x_test,1);

            % remove dimensions with no variance
            dim_oi = var(x_train) > 0;
            x_train = x_train(:, dim_oi);
            x_test = x_test(:, dim_oi);
            x_ob = x_ob(:, dim_oi);

            % standardize dimensions
            mu = nanmean(x_train);
            sig = nanstd(x_train);
            x_ob = (x_ob - mu) ./ sig;
            x_train = (x_train - mu) ./ sig;
            x_test = (x_test - mu) ./ sig;
        
            % rotate to pc axes
            if run_pca 
                [coeff, ~, latent] = pca(x_train);
                varexp = cumsum(latent)./sum(latent);
%                 t_pc = find(varexp >= 1, 1, 'first');
                t_pc = size(x_train,2);
                x_train = x_train * coeff;  x_train = x_train(:,1:t_pc);
                x_test = x_test * coeff;    x_test = x_test(:,1:t_pc);
                x_ob = x_ob * coeff;      x_ob = x_ob(:,1:t_pc);
            end    
        
            % compute pairwise distances
            D = squareform(pdist([x_ob; x_test; x_train], dist_met)).^2;
            tmp = D(1,ntest+2:end);
            d_0(iter,ki,1:length(tmp)) = tmp;
            tmp = D(2:(ntest+1), ntest+2:end);
            d_dist(iter,te_inds, 1:size(tmp,2)) = tmp;
        end
    end
        
    d_0_ = nanmean(d_0(:)); 
    d_dist_ = nanmean(nanmean(d_dist,3),1);
    d_dist_ = d_dist_(:);
    
    pd = determine_distribution(d_dist_, dist_opts);
    del_dist = repmat(d_0_,size(d_dist_,1),1) - d_dist_;
    dist_data.d_dist = d_dist_;
    dist_data.d_0 = d_0_;
    dist_data.del_dist = del_dist;
    p = 1 - cdf(pd, d_0_);
    z = (d_0_ - nanmean(d_dist_)) ./ nanstd(d_dist_);
        
    function pd = determine_distribution(d, dist_opts)
        if ~exist('dist_opts', 'var') || isempty(dist_opts)
            pd = force_kernel_fit(d);
            return;
        end
        xvals = linspace(min(d), max(d), 100);        
        ND = length(dist_opts);
        PDs = cell(ND,1);
        KS_score = nan(ND,1);
        for doi = 1:ND
            PDs{doi} = fitdist(d, dist_opts{doi});
            if ND == 1
                pd = PDs{doi};
                return;
            else
                cdf_val = [xvals(:), cdf(PDs{doi}, xvals(:))];
                [~,~,KS_score(doi)] = kstest(d, 'CDF', cdf_val);
            end
        end
        [~,mi] = min(KS_score);
        pd = PDs{mi};
    end
    
    function pd = force_kernel_fit(d)
        pd = fitdist(d,'Kernel','Kernel','epanechnikov');
    end

end