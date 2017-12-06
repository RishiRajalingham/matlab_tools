function out = rr_NHST_exact_decode(X,Y, opt)
% function p = rr_NHST_exact_decode(x_null, x_obs, mode)
% NULL HYPOTHESIS STATISTIC TESTING for high dimensional data, via decoding
% For arbitrarily high dimensional data, we can infer whether two
% distributions are overlapping by projecting onto the best one-dimensional
% subspace (optimized on held-out samples).

% This is a generalization of exact tests for high dimensional data and 
% inferences. Univariate exact tests methods test the null 
% hypothesis that the two sets of values are the same (overlapping).
% Here, the null hypothesis is the same, except for a particular univariate
% distributions are the same, up to differences that
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
    
    data = cat(1, X, Y);
    label = cat(1, zeros(nsx,1), ones(nsy,1));
    
    % Use some folds to optimize the linear projection. Project held-out
    % fold onto this line to get univariate variable.
    inds = crossvalind('kfold', label, opt.k);
    tr = inds ~= 1;
    dat_tr = data(tr,:);     lbl_tr = label(tr,:);
    dat_te = data(~tr,:);    lbl_te = label(~tr,:);
    
    if strcmp(opt.cls, 'svm')
        Model = fitcsvm(dat_tr,lbl_tr, 'Standardize',true);
        w = (Model.Beta);
        m = dat_te * w;
    elseif strcmp(opt.cls, 'lda')
        Model = fitcdiscr(dat_tr,lbl_tr);
        w = Model.Coeffs(1,2).Linear;
        m = dat_te * w;
    end
        
    X_proj = m(lbl_te == 0);
    Y_proj = m(lbl_te == 1);
    
    % exact test on univariate variables
    del = X_proj - Y_proj;
    
    pd = fitdist(del, opt.dist_type);
    p_ = cdf(pd, 0);
    p = min(p_, 1-p_)*2;
    
    out.p_left = p_;
    out.p_right = 1-p_;
    out.p = p;
    out.X_proj = X_proj;
    out.Y_proj = Y_proj;
    

  

end
