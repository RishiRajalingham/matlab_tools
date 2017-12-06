function out = rr_NHST_decode(X,Y, opt)
% function p = rr_NHST_decode(x_null, x_obs, mode)
% NULL HYPOTHESIS STATISTIC TESTING for high dimensional data, via decoding
% For arbitrarily high dimensional data, we can infer whether two
% distributions are different by attempting to decode the dist label, and
% compare this to shuffled decode using univariate exact test.

% This is a generalization of t-tests, sign-rank tests, etc.., for high
% dimensional data and inferences. Those other methods test the null 
% hypothesis that the means/medians of the two distributions are 
% the same between two sets of observations. Hhere, the null hypothesis is 
% that the generative distributions are the same, up to differences that
% can be observed by the decoders used.

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
