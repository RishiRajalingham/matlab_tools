function [p_,p] = rr_exact_test(x, y, opt)
% function p = rr_exact_test(x, y, opt)
% for univariate variables
    if ~exist('opt', 'var') 
        opt = [];
    end
    
    assert( (size(x,2) == 1) & (size(y,2) == 1));
    del_fn = rr_pop_opts(opt, 'del_fn', 1); % correspondence, or all pairwise?
    density_est = rr_pop_opts(opt, 'density_est', 'kde'); 
    signeddistfun = @(a,b) (a-b);

    if del_fn == 0
        del = signeddistfun(x,y); 
    elseif del_fn == 1
        del = pdist2(x, y, signeddistfun);
    end
    
    del = del(:);
    if strcmp(density_est, 'kde')
        [~,pdf_x,x,~] = kde(del);
        p_ = sum(pdf_x(x<=0))./sum(pdf_x);
    else
        p_ = sum(del < 0)/sum(~isnan(del));
    end
    
    p = min(p_, 1-p_)*2;

end