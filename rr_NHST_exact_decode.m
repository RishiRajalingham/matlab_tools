function out = rr_NHST_exact_decode(X,Y,opt)
% function p = rr_NHST_exact_decode(X, Y, mode)
% NULL HYPOTHESIS STATISTIC TESTING for high dimensional data, via decoding
% For arbitrarily high dimensional data, we can infer whether two
% distributions are overlapping by projecting onto the best one-dimensional
% subspace (optimized on *held-out* data).
%
% This is a generalization of exact tests for high dimensional data and 
% inferences. Univariate exact tests methods test the null 
% hypothesis that the two sets of values are the same (overlapping).
% Here, the null hypothesis is the same, except for a particular univariate
% distributions are the same, up to differences that
% can be observed by the decoders used.
%
% X,Y are cells of size (niter,2), where the second dimension corresponds
% to estimates from split halves of *independent* observations.
% opt has fieldnames: 
%   'cls': ('lda'/'svm') for classifier optimization
%   'del_fn': all pairwise deltas? 0/1

    % Parse options
    if ~exist('opt', 'var') 
        opt = [];
    end
    cls = rr_pop_opts(opt, 'cls', 'lda');
    
    % Ensure that train/test folds exist
    k = 2;
    assert( iscell(X) & iscell(Y) );
    assert( size(X,1) == size(Y,1) );   niter = size(X,1);
    assert( (size(X,2) == k) & (size(Y,2) == k) );
    assert( (size(X{1},2) == size(Y{1},2)) ); nd = size(X{1},2);
    
    % Run
    [X_proj, Y_proj, W] = crossval_linear_remap(X,Y);
    [p_,p] = rr_exact_test(X_proj, Y_proj, opt);
    out = [];
    out.p = p;
    out.p_left = p_;
    out.p_right = 1-p_;
    out.X_proj = X_proj;
    out.Y_proj = Y_proj;
    out.W = W;
    
    function [X_proj, Y_proj, W] = crossval_linear_remap(X,Y)
        X_proj = [];
        Y_proj = [];
        W = nan(niter,k,nd);
        for iter = 1:niter
            for ki = 1:k
                kj = k-ki+1;
                dat_tr = cat(1, X{iter,ki}, Y{iter,ki});
                lbl_tr = cat(1, zeros(size(X{iter,ki},1),1), ones(size(Y{iter,ki},1),1));
                dat_te = cat(1, X{iter,kj}, Y{iter,kj});
                lbl_te = cat(1, zeros(size(X{iter,kj},1),1), ones(size(Y{iter,kj},1),1));
                
                if strcmp(cls, 'svm')
                    Model = fitcsvm(dat_tr,lbl_tr, 'Standardize',true);
                    w = (Model.Beta);
                    m = dat_te * w;
                elseif strcmp(cls, 'lda')
                    Model = fitcdiscr(dat_tr,lbl_tr);
                    w = Model.Coeffs(1,2).Linear;
                    m = dat_te * w;
                end
                
                X_proj = cat(1, X_proj, m(lbl_te == 0));
                Y_proj = cat(1, Y_proj, m(lbl_te == 1));
                W(iter,ki,:) = w;
            end
        end
    end
  
    

end
