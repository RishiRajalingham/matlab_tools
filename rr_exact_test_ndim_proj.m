function p = rr_exact_test_ndim_proj(x_null, x_obs, mode)
% function p = rr_exact_test_ndim(x_null, x_obs)
% For a n-dimensional random variable, get the empirical probability that
% it was sampled from the null distribution, using a 1-d projection method.

% based on simulations, this method is wrong!!!!!
% i.e. it will bias p-values such that false detection rate scales with the
% number of dimension of the variable, i.e. does not correct for multiple 
% comparisons. this is probable because the projection is optimized for the 
% true observation. if done correctly, the null dist should be on random 
% projections. but this is too expensive and not worth it.
    

    
    [ns,nd] = size(x_null);
    if ~exist('x_obs', 'var') || isempty(x_obs)
        x_obs = zeros(1,nd);
    end
     
    % split into disjoint data used for mapping and testing.
    x_train = x_null(1:floor(ns/2),:);
    x_test = x_null(1+floor(ns/2):end,:);
    
    % recenter to mean of null distribution
    x_mu = nanmean(x_train);
%     x_train = x_train - x_mu;
    x_test = x_test - x_mu;
    x_obs = x_obs - x_mu;
    
    % map onto one axis
    x_obs_n = repmat(x_obs,size(x_test,1),1);
    null_mapped = dot(x_test',x_obs_n')' ./ dot(x_obs_n', x_obs_n')';
    obs_mapped = dot(x_obs_n',x_obs_n')' ./ dot(x_obs_n', x_obs_n')';
    del_mapped = null_mapped - obs_mapped; %  obs is at 1
   
    p = rr_normal_exact_test(del_mapped, [], mode);


end
