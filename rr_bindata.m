function [X,Y] = rr_bindata(x, n, nx)
% function [x2,y] = rr_bindata(x, n)

    if ~isempty(n) && isscalar(n)
        xnan = x(~isnan(x));
        x2 = prctile(xnan,0:n:100);
    elseif ~isempty(n) && isvector(n)
        x2 = prctile(x,n);
    elseif ~isempty(nx)
        x2 = nx;
    end
    X = nan(size(x));
    Y = nan(size(x));
    for i = 1:length(x2)-1
        t = x >= x2(i) & x <= x2(i+1);
        Y(t) = i;
        X(t) = mean(x2(i:i+1));
    end
%     [~,t] = max(x);
%     Y(t) = i;
%     X(t) = mean(x2(i:i+1));
end