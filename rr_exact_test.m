function [p_,p] = rr_exact_test(xx, yy, opt)
% function p = rr_exact_test(x, y, opt)
% for univariate variables
% p_ =  (xx-yy) < 0

    if ~exist('opt', 'var') 
        opt = [];
    end
    
    nd = size(xx,2);
    p_ = nan(nd,1);
    p = nan(nd,1);
    for ndi = 1:nd
        if isnan(nanmean(xx(:,ndi))) || isnan(nanmean(yy(:,ndi)))
            continue;
        end
        [p_(ndi),p(ndi)] = exact_test_univariate(xx(:,ndi),yy(:,ndi));
    end
    
    
    function [p_,p] = exact_test_univariate(x,y)
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
            mm = max(abs(del))*1.5;
            [~,pdf_x,x_ax,~] = kde(del, 2^14, -mm, mm);
            p_ = sum(pdf_x(x_ax<=0))./sum(pdf_x);
        else
            p_ = sum(del < 0)/sum(~isnan(del));
        end
%         p_ = p_; %one-tailed
        p = min(p_, 1-p_)*2;
    end

end