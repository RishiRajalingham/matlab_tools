function metric = behaviouralMetrics(data, func)
%  function metric = behaviouralMetrics(data, func)- rishir
% 0 : performance
% 1 : cmatvec
% 2 : dprime
% 3 : asymmetry
% 4 : dprimemat
% 5 : dprimeova
% 6 : semantic
% 7 : cmat
% 8 : bias
% 9 : perfova
% 10 : rt
% 11 : percorr
% 12 : pid


%% Route image-level metrics here
    if ~isempty(strfind(func, 'imglvl'))
        metric = imglvl_behaviouralMetrics(data, func);
        return;
    end

%% Init
    s = data(:,1); % sample category
    m = data(:,2); % match test category
    nm = data(:,3); % nonmatch test category
    sel = data(:,4); % selected category/similarity rating
    if size(data,2) > 4;
        rt = data(:,5); % reaction time/sometimes imgid
    end
    
    if false
        us = unique(s); % all categories
        Ns = length(us); % number of categories
    else
        us = 0:24;
        Ns = 25;
    end
   
%% Parse
    if ischar(func)
        func = getBehMetric_code(func);
    end

    % dprime ceiling
    maxVal = 5;
%% Compute
    switch func
        case getBehMetric_code('performance')
            metric = sum(s == sel) / length(s); 
            return;
        case getBehMetric_code('dprime_mu')
            metric = behaviouralMetrics(data, 'dprime');
            metric = nanmean(metric);
            return;
        case getBehMetric_code('bias_mu')
            metric = behaviouralMetrics(data, 'bias');
            metric = nanmean(abs(metric));
            return;
        case getBehMetric_code('cmatvec')
            di = logical(eye(Ns)); 
            C1 = rr_confusionmat(s,sel, us,us);
            C2 = rr_confusionmat(s,m, us,us) + rr_confusionmat(s,nm, us,us);
            C = C1 ./ C2;
            metric = C(~di);  
            metric = 1 - metric;  
            return;
        case getBehMetric_code('cmatvec_logit')
            di = logical(eye(Ns)); 
            C1 = rr_confusionmat(s,sel, us,us);
            C2 = rr_confusionmat(s,m, us,us) + rr_confusionmat(s,nm, us,us);
            C = C1 ./ C2;
            metric = 1 - C(~di);  
            t = isnan(metric);
            metric = log(metric./(1-metric));
            metric = min(metric, 10);
            metric(t) = nan;
            return;
        case getBehMetric_code('dprime')
            dPrime = ones(Ns*(Ns-1)/2,1).*NaN;
            lapConstant = 0;%1^-10; %laplace smoothing?
            
            count = 0;
            for i = 1:Ns
                si = us(i);
                for j = i+1:Ns
                    sj = us(j);
                    tr = (m == si & nm== sj) | (m == sj & nm == si);
                    count = count+1;
                    if sum(tr) == 0
                        dPrime(count) = NaN;
                    else
                        hr = (sum(sel(tr) == si & s(tr) == si) + lapConstant) /...
                            (sum(s(tr) == si) + 2*lapConstant);
                        fp = (sum(sel(tr) == si & s(tr) == sj) + lapConstant) /...
                            (sum(s(tr) == sj) + 2*lapConstant);
                        dp = norminv(hr,0,1) - norminv(fp,0,1);
                        dPrime(count) = constrain_range(dp, -maxVal, maxVal);
                    end
                end
            end
            metric = dPrime; 
            return;
            
        case getBehMetric_code('balacc')
            balacc = ones(Ns*(Ns-1)/2,1).*NaN;
            lapConstant = 0;%1^-10; %laplace smoothing?
            
            count = 0;
            for i = 1:Ns
                si = us(i);
                for j = i+1:Ns
                    sj = us(j);
                    tr = (m == si & nm== sj) | (m == sj & nm == si);
                    count = count+1;
                    if sum(tr) == 0
                        balacc(count) = NaN;
                    else
                        hr = (sum(sel(tr) == si & s(tr) == si) + lapConstant) /...
                            (sum(s(tr) == si) + 2*lapConstant);
                        fp = (sum(sel(tr) == si & s(tr) == sj) + lapConstant) /...
                            (sum(s(tr) == sj) + 2*lapConstant);
                        balacc(count)  = hr - fp;
                    end
                end
            end
            metric = balacc; 
            return;
        case getBehMetric_code('asymmetry')        
            C1 = confusionmat(s,sel);
            C2 = confusionmat(s,m) + confusionmat(s,nm);
            C = C1 ./ C2;
            metric = norm(C-C');
            return
        case getBehMetric_code('dprimemat')        
            dPrimeMat = nan(Ns,Ns);
%             maxVal = 5;
            lapConstant = 0;

            for i = 1:Ns
                si = us(i);
                for j = i+1:Ns
                    sj = us(j);
                    tr = (m == si & nm== sj) | (m == sj & nm == si);
                    if sum(tr) == 0
                        continue;
                    else
                        hr = (sum(sel(tr) == si & s(tr) == si) + lapConstant) /...
                            (sum(s(tr) == si) + 2*lapConstant);
                        fp = (sum(sel(tr) == si & s(tr) == sj) + lapConstant) /...
                            (sum(s(tr) == sj) + 2*lapConstant);
                        dp = min(norminv(hr,0,1) - norminv(fp,0,1), maxVal);
                        dPrimeMat(i,j) = dp;
                        dPrimeMat(j,i) = dp;
                    end
                end
            end
            metric = dPrimeMat;
            return;
        case getBehMetric_code('dprimeova') % one versus all sensitivity index    
            dPrime = nan(Ns,1);
            lapConstant = 0;%1^-10; %laplace smoothing?
%             maxVal = 5;

            for i = 1:Ns
                si = us(i);
                tr = (m == si) | (nm == si);
                if sum(tr) ~= 0
                    hr = (sum(sel(tr) == si & s(tr) == si) + lapConstant) /...
                        (sum(s(tr) == si) + 2*lapConstant);
                    fp = (sum(sel(tr) == si & s(tr) ~= si) + lapConstant) /...
                        (sum(s(tr) ~= si) + 2*lapConstant);
                    dPrime(i) = min(norminv(hr,0,1) - norminv(fp,0,1), maxVal);
                end

            end
            metric = dPrime;
            return;
        case getBehMetric_code('balaccova') % one versus all sensitivity index    
            dPrime = nan(Ns,1);
            lapConstant = 0;%1^-10; %laplace smoothing?
%             maxVal = 5;

            for i = 1:Ns
                si = us(i);
                tr = (m == si) | (nm == si);
                if sum(tr) ~= 0
                    hr = (sum(sel(tr) == si & s(tr) == si) + lapConstant) /...
                        (sum(s(tr) == si) + 2*lapConstant);
                    fp = (sum(sel(tr) == si & s(tr) ~= si) + lapConstant) /...
                        (sum(s(tr) ~= si) + 2*lapConstant);
                    dPrime(i) = hr - fp;
                end

            end
            metric = dPrime;
            return;
        case getBehMetric_code('biasova') % one versus all bias
            bias = nan(Ns,1);
            lapConstant = 0;%1^-10; %laplace smoothing?
%             maxVal = 5;

            for i = 1:Ns
                si = us(i);
                tr = (m == si) | (nm == si);
                if sum(tr) ~= 0
                    hr = (sum(sel(tr) == si & s(tr) == si) + lapConstant) /...
                        (sum(s(tr) == si) + 2*lapConstant);
                    fp = (sum(sel(tr) == si & s(tr) ~= si) + lapConstant) /...
                        (sum(s(tr) ~= si) + 2*lapConstant);
                    tmp = (norminv(hr,0,1) + norminv(fp,0,1))./2;
                    bias(i) = max(min(tmp, maxVal), -maxVal);
                end

            end
            metric = (bias);
            return;
        case getBehMetric_code('semantic') % semantic similarity rating
            if Ns ~= 24; Ns = 24; end; % overwrite
            semsil = ones(Ns*(Ns-1)/2,1).*NaN;
            count = 0;
            for i = 1:Ns
                si = i-1;%us(i);
                for j = i+1:Ns
                    sj = j-1;%us(j);
                    tr = (m == si & nm== sj) | (m == sj & nm == si);
                    count = count+1;
                    if sum(tr) == 0
                        semsil(count) = NaN;
                    else
                        semsil(count) = mean(sel(tr));
                    end
                end
            end
            metric = semsil;
            return
        case getBehMetric_code('cmat')% pseudo confusion matrix
            di = logical(eye(Ns)); 
            C1 = confusionmat(s,sel);
            C2 = confusionmat(s,m) + confusionmat(s,nm);
            C = C1 ./ C2;
    %         C(di) = nan;
            metric = C;
        case getBehMetric_code('bias')
            bias = ones(Ns*(Ns-1)/2,1).*NaN;
            lapConstant = 0;%1^-10; %laplace smoothing?
            count = 0;
            for i = 1:Ns
                si = us(i);
                for j = i+1:Ns
                    sj = us(j);
                    tr = (m == si & nm== sj) | (m == sj & nm == si);
                    count = count+1;
                    if sum(tr) == 0
                        bias(count) = NaN;
                    else
                        hr = (sum(sel(tr) == si & s(tr) == si) + lapConstant) /...
                            (sum(s(tr) == si) + 2*lapConstant);
                        fp = (sum(sel(tr) == si & s(tr) == sj) + lapConstant) /...
                            (sum(s(tr) == sj) + 2*lapConstant);
                        c = (norminv(hr,0,1) + norminv(fp,0,1)) * 0.5;
                        bias(count) = constrain_range(c, -5, 5);
                    end
                end
            end
            metric = bias;
            return;
            
        case getBehMetric_code('bias2') %asymmetry im hit rates, without non-linearity
            bias = ones(Ns*(Ns-1)/2,1).*NaN;
            lapConstant = 0;%1^-10; %laplace smoothing?
            count = 0;
            for i = 1:Ns
                si = us(i);
                for j = i+1:Ns
                    sj = us(j);
                    tr = (m == si & nm== sj) | (m == sj & nm == si);
                    count = count+1;
                    if sum(tr) == 0
                        bias(count) = NaN;
                    else
                        hr = (sum(sel(tr) == si & s(tr) == si) + lapConstant) /...
                            (sum(s(tr) == si) + 2*lapConstant);
                        fp = (sum(sel(tr) == si & s(tr) == sj) + lapConstant) /...
                            (sum(s(tr) == sj) + 2*lapConstant);
                        bias(count) = mean([hr,fp]);
                    end
                end
            end
            metric = bias;
            return;
        case getBehMetric_code('perfova') % one versus all percent correct
            perf = nan(Ns,1);
            for i = 1:Ns
                si = us(i);
                tr = (m == si) | (nm == si);
                if sum(tr) ~= 0
                    perf(i) = sum(sel(tr) == s(tr)) / sum(tr);
                end

            end
            metric = perf;
            return
        
        case getBehMetric_code('perfova_logit') % one versus all percent correct
            perf = nan(Ns,1);
            for i = 1:Ns
                si = us(i);
                tr = (m == si) | (nm == si);
                if sum(tr) ~= 0
                    perf(i) = sum(sel(tr) == s(tr)) / sum(tr);
                end

            end
            metric = perf;
            t = isnan(metric);
            metric = log(metric./(1-metric));
            metric = min(metric, 10);
            metric(t) = nan;
            return
            
        case getBehMetric_code('rt')% reaction time per binary task
            reacTime = ones(Ns*(Ns-1)/2,1).*NaN;
            maxVal = 15*10^6;

            count = 0;
            for i = 1:Ns
                si = us(i);
                for j = i+1:Ns
                    sj = us(j);
                    tr = (m == si & nm== sj) | (m == sj & nm == si);
                    count = count+1;
                    rt_tmp = rt(tr);
                    reacTime(count) = nanmean(rt_tmp(rt_tmp < maxVal));
                end
            end
            metric = reacTime; 
            return;
        case getBehMetric_code('rtova')% reaction time per binary task
            reacTime = ones(Ns,1).*NaN;
            maxVal = 15*10^6;
            for i = 1:Ns
                si = us(i);
                tr = (m == si) | (nm == si);
                rt_tmp = rt(tr);
                reacTime(i) = nanmean(rt_tmp(rt_tmp < maxVal));
            end
            metric = reacTime; 
            return;
        
        case getBehMetric_code('percorr')% percent correct
            percorr = nan(Ns*(Ns-1)/2,1);
            count = 0;
            for i = 1:Ns
                si = us(i);
                for j = i+1:Ns
                    sj = us(j);
                    tr = (m == si & nm== sj) | (m == sj & nm == si);
                    count = count+1;
                    if sum(tr) ~= 0
                        percorr(count) = (sum(sel(tr) == s(tr))) / (sum(tr)); 
                    end
                end
            end
            metric = percorr;
            return;
            
        case getBehMetric_code('percorr_logit')% percent correct
            percorr = nan(Ns*(Ns-1)/2,1);
            count = 0;
            for i = 1:Ns
                si = us(i);
                for j = i+1:Ns
                    sj = us(j);
                    tr = (m == si & nm== sj) | (m == sj & nm == si);
                    count = count+1;
                    if sum(tr) ~= 0
                        percorr(count) = (sum(sel(tr) == s(tr))) / (sum(tr)); 
                    end
                end
            end
            metric = percorr;
            t = isnan(metric);
            metric = log(metric./(1-metric));
            metric = min(metric, 10);
            metric(t) = nan;
            return;
            
        case getBehMetric_code('pid')
            pid = nan(Ns*(Ns-1)/2,1);
            count = 0;
            for i = 1:Ns
                si = us(i);
                for j = i+1:Ns
                    sj = us(j);
                    count = count+1;
                    pid(count) = si*100+sj;
                end
            end
            metric = pid;
            return;
        case getBehMetric_code('pid2')
            pid = nan(Ns*(Ns-1),1);
            count = 0;
            for i = 1:Ns
                si = us(i);
                for j = 1:Ns
                    sj = us(j);
                    if si == sj; continue; end;

                    count = count+1;
                    pid(count) = sj*100+si;
                end
            end
            metric = pid;
            return;
    end
    
    
    
    
    
    
    %% Helper 
    function CC = rr_confusionmat(G1, G2, o1, o2)
        [C_,o] = confusionmat(G1, G2);
        t1 = ismember(o,o1);
        t2 = ismember(o,o2);
        s1 = ismember(o1,o);
        s2 = ismember(o2,o);
        CC = nan(length(o1), length(o2));
        CC(s1,s2) = C_(t1,t2);
    end
  
end