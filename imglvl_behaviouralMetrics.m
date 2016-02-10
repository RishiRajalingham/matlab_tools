function metric = imglvl_behaviouralMetrics(data, func)
% function metric = imglvl_behaviouralMetrics(data, func)

%% Init
    s = data(:,1); % sample category
    m = data(:,2); % match test category
    nm = data(:,3); % nonmatch test category
    sel = data(:,4); % selected category/similarity rating
    img = data(:,5); % img-index
    
    us = unique(s); % all categories
    Ns = length(us); % number of categories
    uimg = 1:2400;
    nImg = numel(uimg);
    nimg_perobj = 100;
    metric = [];
       
%% Parse 
    if strcmp(func, 'imglvl_performance') == 1
        func = 0 ;
    elseif strcmp(func, 'imglvl_dprimeova') == 1
        func = 1;
    elseif strcmp(func, 'imglvl_cmat') == 1 
        func = 2;
    elseif strcmp(func, 'imglvl_cmatvec') == 1 
        func = 3;
    elseif strcmp(func, 'imglvl_ntrials') == 1 
        func = 4;
    elseif strcmp(func, 'imglvl_pid') == 1 
        func = 5;
    elseif strcmp(func, 'imglvl_performance_norm') == 1
        func = 6;
    elseif strcmp(func, 'imglvl_dprimeovo') == 1
        func = 7;
    end
   
%% Compute
    if func == 0% image-level performance
        
        succ = s == sel;
        t = ~isnan(img);
        [a,c] = grpstats(succ(t), img(t), {'mean','gname'});
        c = str2double(c);
        mu = nan(nImg,1);
        mu(c) = a;
        metric = mu;
        
    elseif func == 1% image-level unbiased performance (one versus all dp)
        
        C1 = rr_confusionmat(img, sel, uimg, us);
        C2 = rr_confusionmat(img, m, uimg, us) ...
            + rr_confusionmat(img, nm, uimg, us);
        C = C1;
        C(C2 == 0) = nan;
        
        true_obj = floor((uimg-1)./nimg_perobj);
        C3 = rr_confusionmat(uimg, true_obj, uimg, us);
        
        C_rowsum = nansum(C,2);
        
        maxVal = 5;
        unbiased_perf = nan(nImg,1);
        for i = 1:nImg
            cat_i = logical(C3(i,:));
            hr = C(i,cat_i) ./ C_rowsum(i);
            tmp = C3(:,cat_i) == 0;
		if sum(cat_i) == 0 | sum(tmp) == 0
			continue;
		end
            fp = nanmean(C(tmp,cat_i) ./ C_rowsum(tmp));
            up_ = norminv(hr,0,1) - norminv(fp,0,1);
            if ~isnan(up_)
                unbiased_perf(i) = min(up_, maxVal);
            end
        end
        
        metric = unbiased_perf;
        
    elseif func == 2% cmat image-level confusion

        C1 = rr_confusionmat(img, sel, uimg, us);
        C2 = rr_confusionmat(img,m,uimg, us) ...
            + rr_confusionmat(img,nm,uimg, us);
        C = C1./C2;
        C(C2 == 0) = nan;
        
        true_obj = floor((uimg-1)./nimg_perobj);
        C3 = rr_confusionmat(uimg,true_obj,uimg, us);
        C(C3 > 0) = nan;
        metric = C;
    
    elseif func == 3% cmat image-level confusion

        C = imglvl_behaviouralMetrics(data, 'imglvl_cmat');
        metric = C(:);
        
    elseif func == 4% raw cmat image-level confusion

        C = rr_confusionmat(img,nm,uimg, us);
        metric = C;
    elseif func == 5% unique id for each image-cat pair
        pid = img.*1000 + nm;
        metric = pid;
    elseif func == 6% image-performance normalized by object performance
        mu = imglvl_behaviouralMetrics(data, 'imglvl_performance');
        mu_n = mu;
        for ni = 1:Ns
            imgi = (ni-1)*nimg_perobj + 1 : (ni)*nimg_perobj;
            mu_n(imgi) = (mu(imgi) - nanmean(mu(imgi))) ./ nanstd(mu(imgi));
        end
        metric = mu_n;
    elseif func == 7
        C1 = rr_confusionmat(img, sel, uimg, us);
        C2 = rr_confusionmat(img, m, uimg, us) ...
            + rr_confusionmat(img, nm, uimg, us);
        C = C1;
        C(C2 == 0) = nan;
        
        true_obj = floor((uimg-1)./nimg_perobj);
        C3 = rr_confusionmat(uimg, true_obj, uimg, us);
       
        unbiased_perf = nan(nImg, Ns);
        for i = 1:nImg
            % which category is this image
            cat_i = find(C3(i,:) == 1); 
            for nsi = 1:Ns
                if nsi == cat_i; continue; end;
                % what images are not in this category
                ncat_i = C3(:,cat_i) == 0; 
                C2x2 = [C(i,[cat_i,nsi]); C(ncat_i,[cat_i, nsi])];
                unbiased_perf(i,nsi) = dprime_from2x2(C2x2); 
            end
        end
        
        metric = unbiased_perf(:);
    end
    
    
    %% Helper functions
    function dp = dprime_from2x2(C2x2)
        maxVal = 5;
        hr_ = C2x2(1,1) ./ nansum(C2x2(1,:),2);
        fp_ = nansum(C2x2(2:end,1)) ./ nansum(nansum(C2x2(2:end,:),2));
        dp = norminv(hr_,0,1) - norminv(fp_,0,1);
        dp = min(dp, maxVal);
    end
    
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
