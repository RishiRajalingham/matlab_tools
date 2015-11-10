function N = getNumSplitHalves(N1)
% function getCorrelationType(ctype)
% no arguments for get, 1 argument for set

    persistent RR_NSPLITHALVES
    if nargin == 1
        RR_NSPLITHALVES = N1;
    else
        N = RR_NSPLITHALVES;
    end
end