function ctype = getCorrelationType(ctype1)
% function getCorrelationType(ctype)
% no arguments for get, 1 argument for set

    persistent RR_CORRELATION_TYPE
    if nargin == 1
        RR_CORRELATION_TYPE = ctype1;
    else
        ctype = RR_CORRELATION_TYPE;
    end
end