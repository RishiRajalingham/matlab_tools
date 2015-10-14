function metric = behaviouralMetrics(data, func)
%  function metric = behaviouralMetrics(data, func)
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

    s = data(:,1); % sample category
    m = data(:,2); % match test category
    nm = data(:,3); % nonmatch test category
    sel = data(:,4); % selected category/similarity rating
    
    if size(data,2) > 4;
        rt_img = data(:,5); % reaction time or img-index
    end
    
    us = unique(s); % all categories
    Ns = length(us); % number of categories
    
    if strcmp(func, 'performance') == 1
        func = 0 ;
    elseif strcmp(func, 'cmatvec') == 1
        func = 1;
    elseif strcmp(func, 'dprime') == 1 
        func = 2;
    elseif strcmp(func, 'asymmetry') == 1
        func = 3;
    elseif strcmp(func, 'dprimemat') == 1 
        func = 4;
    elseif strcmp(func, 'dprimeova') == 1 
        func = 5;
    elseif strcmp(func, 'semantic') == 1
        func = 6;
    elseif strcmp(func, 'cmat') == 1
        func = 7;
    elseif strcmp(func, 'bias') == 1
        func = 8;
    elseif strcmp(func, 'perfova') == 1
        func = 9 ;
    elseif strcmp(func, 'rt') == 1
        func = 10;
    elseif strcmp(func, 'percorr') == 1
        func = 11;
    elseif strcmp(func, 'pid') == 1
        func = 12;
    elseif strcmp(func, 'pid2') == 1
        func = 13;
    elseif strcmp(func, 'imgperf') == 1
        func = 14;
    elseif strcmp(func, 'imgdprime') == 1
        func = 15;
    elseif strcmp(func, 'imgconf') == 1
        func = 16;
    end
    
    if func == 0 % performance
        metric = sum(s == sel) / length(s);
        
    elseif func == 0.5 % percent correct over tasks
        pcorr = ones(Ns*(Ns-1)/2,1).*NaN;
        
        count = 0;
        for i = 1:Ns
            si = us(i);
            for j = i+1:Ns
                sj = us(j);
                tr = (m == si & nm== sj) | (m == sj & nm == si);
                count = count+1;
                if sum(tr) == 0
                    pcorr(count) = NaN;
                else
                    pcorr(count) = sum(sel(tr) == s(tr)) / sum(tr);
                end
            end
        end
        metric = pcorr;
        
    elseif func == 1% pseudo confusion matrix - unwrapped
        
        di = logical(eye(Ns)); 
        C1 = confusionmat(s,sel);
        C2 = confusionmat(s,m) + confusionmat(s,nm);
        C = C1 ./ C2;
        metric = C(~di);  
        
    elseif func == 2% sensitivity index
        
        dPrime = ones(Ns*(Ns-1)/2,1).*NaN;
        lapConstant = 0;%1^-10; %laplace smoothing?
        maxVal = 5;
        
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
                    dPrime(count) = min(norminv(hr,0,1) - norminv(fp,0,1), maxVal);
                end
            end
        end
        metric = dPrime;
        
    elseif func == 3% asymmetry index
        
        C1 = confusionmat(s,sel);
        C2 = confusionmat(s,m) + confusionmat(s,nm);
        C = C1 ./ C2;
        metric = norm(C-C');
    
    elseif func == 4 % dprime matrix
        dPrimeMat = nan(Ns,Ns);
        maxVal = 5;
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
        
    elseif func == 5% one versus all sensitivity index
        
        dPrime = nan(Ns,1);
        lapConstant = 0;%1^-10; %laplace smoothing?
        maxVal = 5;
        
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
        
    elseif func == 6% semantic similarity rating
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
       
    elseif func == 7 % pseudo confusion matrix
        di = logical(eye(Ns)); 
        C1 = confusionmat(s,sel);
        C2 = confusionmat(s,m) + confusionmat(s,nm);
        C = C1 ./ C2;
%         C(di) = nan;
        metric = C;
        
    elseif func == 8 % bias for each task
        
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
                    bias(count) = max(min((norminv(hr,0,1) + norminv(fp,0,1)) * 0.5, 5), -5);
                end
            end
        end
        metric = bias;
    elseif func == 9 % one versus all percent correct
        
        perf = nan(Ns,1);
        for i = 1:Ns
            si = us(i);
            tr = (m == si) | (nm == si);
            if sum(tr) ~= 0
                perf(i) = sum(sel(tr) == s(tr)) / sum(tr);
            end
        
        end
        metric = perf;
    elseif func == 10 % reaction time per binary task
        
        reacTime = ones(Ns*(Ns-1)/2,1).*NaN;
        maxVal = 15*10^6;
        
        count = 0;
        for i = 1:Ns
            si = us(i);
            for j = i+1:Ns
                sj = us(j);
                tr = (m == si & nm== sj) | (m == sj & nm == si);
                count = count+1;
                rt_tmp = rt_img(tr);
                reacTime(count) = nanmean(rt_tmp(rt_tmp < maxVal));
            end
        end
        metric = reacTime;  
    elseif func == 11% percent correct
        
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
    elseif func == 12% pair id
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
    
    elseif func == 13% pair id asym
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
        
    elseif func == 14% image-level performance
        
        img = data(:,5);
        succ = data(:,1) == data(:,4);
        t = ~isnan(img);
        [a,c] = grpstats(succ(t), img(t), {'mean','gname'});
        c = str2double(c);
        mu = nan(2400,1);
        mu(c) = a;
        metric = mu;
        
    elseif func == 15% image-level unbiased performance
        uimg = 1:2400;
        uobj = 0:Ns-1;
        nImg = numel(uimg);
        nimg_percat = 100;
        img_i = rt_img;
        
        C1 = rr_confusionmat(img_i, sel, uimg, uobj);
        C2 = rr_confusionmat(img_i,m,uimg, uobj) ...
            + rr_confusionmat(img_i,nm,uimg, uobj);
        C = C1;
        C(C2 == 0) = nan;
        
        true_obj = floor((uimg-1)./nimg_percat);
        C3 = rr_confusionmat(uimg,true_obj,uimg, uobj);
        
        C_rowsum = nansum(C,2);
        
        maxVal = 10;
        unbiased_perf = nan(nImg,1);
        for i = 1:nImg
            cat_i = logical(C3(i,:));
            hr = C(i,cat_i) ./ C_rowsum(i);
            tmp = C3(:,cat_i) == 0;
            fp = nanmean(C(tmp,cat_i) ./ C_rowsum(tmp));
            unbiased_perf(i) = min(norminv(hr,0,1) - norminv(fp,0,1), maxVal);
        end
        
        metric = unbiased_perf;
        
    elseif func == 16% image-level confusion
        uimg = 1:2400;
        uobj = 0:Ns-1;
        nImg = numel(uimg);
        nimg_percat = 100;
        img_i = rt_img;
        
        C1 = rr_confusionmat(img_i, sel, uimg, uobj);
        C2 = rr_confusionmat(img_i,m,uimg, uobj) ...
            + rr_confusionmat(img_i,nm,uimg, uobj);
        C = C1./C2;
        C(C2 == 0) = nan;
        
        true_obj = floor((uimg-1)./nimg_percat);
        C3 = rr_confusionmat(uimg,true_obj,uimg, uobj);
        
        selector = true(size(C));
        selector(C3 > 0) = false;
        
        C = C(selector);
        metric = C(:);
        
    end
    
    
    %% Helper functions
    function C = rr_confusionmat(G1, G2, o1, o2)
        [C_,o] = confusionmat(G1, G2);
        t1 = find(ismember(o,o1));
        t2 = find(ismember(o,o2));
        C = C_(t1,t2);
    end

end