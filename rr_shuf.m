function x2 = rr_shuf(x)
    x2 = x(randperm(length(x)),:);
end