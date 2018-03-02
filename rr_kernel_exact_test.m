function [p,p1,p2] = rr_kernel_exact_test(x_null, x_sample)
% function [p,p1,p2] = rr_kernel_exact_test(x_null, x_sample)

    range_val = max(abs([x_null(:); x_sample])).*10;
    np = 10^4;
    xrange = linspace(-range_val, range_val, np);

    [f,xi] = ksdensity(x_null, xrange);
    p1 = nansum(f(xi >= x_sample)) ./ nansum(f);
    p2 = 1-p1;
    p = min(p1, p2) .* 2;
    
   

end