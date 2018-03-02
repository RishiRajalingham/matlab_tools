function [hl, hp] = rr_boundedline(x,y,e, col)
% h = function rr_boundedline(x,y,e)
% simplified wrapped for boundedline which automatically calculates the
% errorbar variable from the SD/SE variable.

        nt = length(x);
        if size(y,2) ~= nt && size(y,1) == nt
            y = y';
            e = e';
        end
        
        if ~exist('col', 'var')
            col = 'b';
        end
        
        sz = size(y);
        nl = sz(1); 
        
        errb = nan(nt,2,nl);
        errb(:,1,:) = e'; %./2;
        errb(:,2,:) = e'; %./2;
        
        
        [hl, hp] = boundedline(x, y, errb, col);
        

end