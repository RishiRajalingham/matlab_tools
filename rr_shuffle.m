function x_s = rr_shuffle(x, ctrvar)

    if ~exist('ctrvar', 'var')
        ctrvar = ones(size(x));
    end

    x_s = x;
    uc = unique(ctrvar);
    for uci = 1:length(uc)
        ts = find(ctrvar == uc(uci));
        rp = randperm(length(ts));
        ts2 = ts(rp);
        x_s(ts,:) = x(ts2,:);
    end

end