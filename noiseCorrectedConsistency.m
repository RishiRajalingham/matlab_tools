function stats = noiseCorrectedConsistency(DATA1, DATA2, metric, terse, STATS1, STATS2)
    
    corrtype = getCorrelationType();
    if isempty(corrtype); corrtype = 'spearman'; end;
    stats = [];
    
    [a3, m1, m2] = raw_correlation(DATA1, DATA2, metric);
    stats.rho.mu = a3;
    stats.rho.sig = [];
    stats.rho.dist = [];
    
    stats.metric.dat1 = m1;
    stats.metric.dat2 = m2;
    
    if ~terse % recompute IC
        [a1,b1,c1] = splitHalfInternalConsistency(DATA1, metric);
        stats.IC1.mu = a1;
        stats.IC1.sig = b1;
        stats.IC1.dist = c1;
        [a2,b2,c2] = splitHalfInternalConsistency(DATA2, metric);
        stats.IC2.mu = a2;
        stats.IC2.sig = b2;
        stats.IC2.dist = c2;
    else % use precomputed internal consistency distributions
        stats.IC1 = STATS1.IC;
        stats.IC2 = STATS2.IC;
    end
    
    try
        IC1 = randsample(stats.IC1.dist, 100, 'true');
        IC2 = randsample(stats.IC2.dist, 100, 'true');
        hc_dist = real(stats.rho.mu ./ sqrt(IC1.*IC2));
        hc_dist = hc_dist(isfinite(hc_dist) & (hc_dist ~= 0));
        stats.rho_n.dist = hc_dist;
        stats.rho_n.mu = nanmean(hc_dist);
        stats.rho_n.sig = nanstd(hc_dist);
    catch e
    end
    
    function [rho, m1, m2] = raw_correlation(data1, data2, metricn)
        metric_1 = behaviouralMetrics(data1, metricn); 
        metric_2 = behaviouralMetrics(data2, metricn); 
        t = isfinite(metric_1) & isfinite(metric_2);
        rho = corr(metric_1(t), metric_2(t), 'type', corrtype);
        m1 = metric_1;
        m2 = metric_2;
    end
 

end