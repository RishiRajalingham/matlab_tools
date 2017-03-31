function [p,z] = rr_normal_exact_test(x_null, x_obs, mode)
% function [p,z] = rr_normal_exact_test(x_null, x_obs, mode)
% mode: (1: one tailed negative), (2: one tailed positive), (0: two tailed)
% for normally distributed variables (e.g. samples of i.i.d), fit gaussian
% to distribution and estimate p value from this functional form.

    if ~exist('mode', 'var')
        mode = 0;
    end
    if ~exist('x_obs', 'var') || isempty(x_obs)
        x_obs = zeros(size(x_null));
    end
    if norm(size(x_obs) - size(x_null)) ~= 0
        del = repmat(x_obs, size(x_null,1),1) - x_null;
    else
        del = x_obs - x_null;
    end
 
    p2 = normcdf(0, nanmean(del), nanstd(del)); % prob of 0 less than mean(x)
    p1 = 1-p2; % prob of greater than mean(x)
    
    z = nanmean(del)./nanstd(del); % how many SDs away from zero?
    
    if mode == 1 % one-tailed negative
        p = p1;
    elseif mode == 2 % one tailed positive
        p = p2;
    elseif mode == 0 % two tailed
        p = min(p1,p2) *2; 
    end
    
    t0 = nanstd(del) == 0;
    p(t0) = nan;
    
   

end