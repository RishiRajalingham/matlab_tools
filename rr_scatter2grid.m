function [Z,X,Y] = rr_scatter2grid(xy,z, smoothfactor)

    if nargin < 3
        smoothfactor = [];
    end

    offset = 1-min(floor(xy(:)));
    X = floor(xy(:,1))+offset;
    Y = floor(xy(:,2))+offset;
    
    Z = nan(max(X),max(X));
    ind = sub2ind(size(Z),X,Y);
    ind = unique(ind);
    Z(ind) = z;
    Z = flipud(Z);
    
end

