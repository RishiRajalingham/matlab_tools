function [w,err] = rr_lsqnonneg_lasso(X,y,lambda)
% function w = rr_lsqnonneg_lasso(X,y,lambda)
% non-negative least squares with lasso regularization 
% https://scicomp.stackexchange.com/questions/10671/tikhonov-regularization-in-the-non-negative-least-square-nnls-pythonscipy

if ~exist('lambda', 'var')
    nlambda = 1000;
    lambda_range = linspace(eps,1,nlambda);
    err = nan(nlambda,1);
    for li = 1:length(lambda_range)
        lambda = lambda_range(li);
        [w,err(li)] = nonneg_lasso(X,y, lambda);
        if norm(w) < eps
            break;
        end
    end
    mi = find(abs(err - min(err)) > eps, 1, 'first');
    if isempty(mi)
        mi = 1;
    end
    lambda = lambda_range(mi);
    [w,err] = nonneg_lasso(X,y, lambda);
else
    [w,err] = nonneg_lasso(X,y, lambda);
end
    
    function [w,err] = nonneg_lasso(X,y, lambda)
        [~, ndim] = size(X);
        lI = lambda.*eye(ndim);
        X2 = cat(1, X, lI);
        y2 = cat(1, y, zeros(ndim,1));
        w = lsqnonneg(X2, y2);
        err = norm(y2 - X2*w);
    end
        
end