function out = rr_exact_test_decode(X,Y, opt)
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

    out = [];
    
    if ~exist('Y', 'var') || isempty(Y)
        fprintf(1, 'Error: Y is empty')    
        return;
    end
    [nsx,ndx] = size(X);
    [nsy,ndy] = size(Y);
    if ndx ~= ndy
        fprintf(1, 'Error: X,Y dimension mismatch')    
        return;
    end
    
    if ~exist('opt', 'var')
        opt.cls = 'lda';
        opt.k = 2;
        opt.nshuf = 100;
        opt.dist_type = 'Normal';
    end

     
    %% Format data
    data = cat(1, X, Y);
    label = cat(1, zeros(nsx,1), ones(nsy,1));

    %% Decode (true/shuffled)
    err = decode(data, label);
    err_s = nan(opt.nshuf,1);
    for si = 1:opt.nshuf
        label_shuf = label(randperm(length(label)));
        err_s(si) = decode(data, label_shuf);
    end
    
    %% Convert to p-value
    pd = fitdist(err_s, opt.dist_type);
    p = cdf(pd, err);
    
    out.p = p;
    out.err = err;
    out.err_s = err_s;
    
    %% Helpers
    function loss = decode(dat, lbl)
        if strcmp(opt.cls, 'svm')
            loss = decode_svm(dat, lbl);
        elseif strcmp(opt.cls, 'lda')
            loss = decode_lda(dat, lbl);
        end
    end
    
    function loss = decode_svm(dat, lbl)
        SVMModel = fitcsvm(dat,lbl,'Standardize',false);
        CVSVMModel = crossval(SVMModel, 'KFold', opt.k);
        loss = kfoldLoss(CVSVMModel);
    end

    function loss = decode_lda(dat, lbl)
        Model = fitcdiscr(dat,lbl);
        CVModel = crossval(Model, 'KFold', opt.k);
        loss = kfoldLoss(CVModel);
    end

end
