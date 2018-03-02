function CC = rr_confusionmat(G1, G2, o1, o2)
% function CC = rr_confusionmat(G1, G2, o1, o2)

    if ~exist('o1', 'var') || isempty(o1)
        o1 = unique(G1);
    end

    if ~exist('o2', 'var') || isempty(o2)
        o2 = unique(G2);
    end
    
    [C_,o] = confusionmat(G1, G2);
    t1 = ismember(o,o1);
    t2 = ismember(o,o2);
    s1 = ismember(o1,o);
    s2 = ismember(o2,o);
    CC = nan(length(o1), length(o2));
    CC(s1,s2) = C_(t1,t2);
end