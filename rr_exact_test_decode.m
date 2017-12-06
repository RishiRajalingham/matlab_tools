function p = rr_exact_test_decode(X,Y, opt)
% function p = rr_exact_test_proj(x_null, x_obs, mode)
% For arbitrarily high dimensional data, we can infer whether two
% distributions are different by attempting to decode the dist label, and
% compare this to shuffled decode using univariate exact test.

% this is appropriate when we have two distributions (e.g. optained by
% bootstrap resampling). not when only one distribution is compared to a
% sample (e.g. permutation tests with shuffled distributions).

% caveat: traditional statistical tests are against the null hypothesis
% that the means/medians of the two distributions are different. here, it's
% not just about the mean, it's about the distribution as a whole. even
% differences in variance will be detected. i.e. these samples were
% generated from different distributions.

    if ~exist('Y', 'var') || isempty(Y)
        fprintf(1, 'Error: Y is empty')
        p = nan;
        return;
    end

    if ~exist('opt', 'var')
        fprintf(1, 'Default opt')
        opt.k = 2;
        opt.nshuf = 100;
        opt.dist_type = 'Normal';
    end

    [nsx,ndx] = size(X);
    [nsy,ndy] = size(Y);
    assert(ndx == ndy);
    k = opt.k;
    nshuf = opt.nshuf;

    data = cat(1, X, Y);
    label = cat(1, zeros(nsx,1), ones(nsy,1));

    % True error rate
    SVMModel = fitcsvm(data,label,'Standardize',true);
    CVSVMModel = crossval(SVMModel, 'KFold', k);
    err = kfoldLoss(CVSVMModel);
    
    % Shuffled error distribution
    err_s = nan(nshuf,1);
    for si = 1:nshuf
        label_shuf = label(randperm(length(label)));
        SVMModel = fitcsvm(data,label_shuf,'Standardize',true);
        CVSVMModel = crossval(SVMModel, 'KFold', k);
        err_s(si) = kfoldLoss(CVSVMModel);
    end
    
    % Convert to p-value
    pd = fitdist(err_s,opt.dist_type);
    p = cdf(pd, err);

end
