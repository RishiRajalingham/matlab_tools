function r = rr_nancorr(x,y,t)
    nn = ~isnan(x) & ~isnan(y);
    if ~exist('t', 'var')
        t = 'spearman';
    end
    if isempty(x(nn)) | isempty(y(nn))
        r = nan;
    else
        r = corr(x(nn), y(nn), 'type', t);
    end
end