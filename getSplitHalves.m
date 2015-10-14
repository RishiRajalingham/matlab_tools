function split_inds = getSplitHalves(data, N, M, imglvl)
% function ind = getSplitHalves(data, N, M)
% get N sets of two balanced splits, 
% number of trials per condition defined by M
    
    m = data(:,2); % match test category
    nm = data(:,3); % nonmatch test category
    imgi = data(:,5); % image index
    
    Ns = length(unique(m)); % number of categories
    Ni = length(unique(imgi)); % number of images
    N_mult = 10^max(ceil(log10(Ns)),ceil(log10(Ni))); % concat multiplier 
    
    obj_task = m.*N_mult + nm;
    img_task = imgi.*N_mult + nm;
    
    if ~exist('imglvl', 'var') || isempty(imglvl) || (imglvl == 0)
        splitby = obj_task;
    elseif imglvl == 1;
        splitby = imgi;
    elseif imglvl == 2;
        splitby = img_task;
    end
    
    u_splitby = unique(splitby);
    n_splitby = length(u_splitby);
    split_inds = cell(N,1);
    for ni = 1:N
        split_inds{ni} = cell(2,1);            
        split_inds{ni}{1} = [];
        split_inds{ni}{2} = [];
    end
       
    for ui = 1:n_splitby
        inds = find(splitby == u_splitby(ui));
        inds = inds(1:min(length(inds), M));
        ninds = length(inds);
        if ninds < 2; continue; end;
        for ni = 1:N
            
            rp = randperm(ninds);
            rp1 = rp(1:floor(ninds/2));
            rp2 = rp(1+floor(ninds/2):end);
            
            split_inds{ni}{1} = vertcat(split_inds{ni}{1}, inds(rp1));
            split_inds{ni}{2} = vertcat(split_inds{ni}{2}, inds(rp2));
        end
    end
    
%     
%         ind = cell(N,1);
%     for ni = 1:N
%         ind{ni} = cell(2,1);
%         ind{ni}{1} = []; ind{ni}{2} = [];
% 
%         for ui1 = 1:Ns
%             for ui2 = 1:Ns
%                 i = find(m == us(ui1) & nm == us(ui2));
%                 i = i(randperm(length(i)));
%                 if exist('M', 'var'); 
%                     i = i(1:min(length(i), M));
%                 end
%                 ind{ni}{1} = [ind{ni}{1}; i(1:floor(length(i)/2))];
%                 ind{ni}{2} = [ind{ni}{2}; i(floor(length(i)/2)+1:end)];
%             end
%         end 
%     end
%     
%     
    
end