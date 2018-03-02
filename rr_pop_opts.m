function val = rr_pop_opts(opts, valname, default_val)
    if isempty(opts)
        val = default_val;
    elseif ~ismember(valname, fieldnames(opts))
        val = default_val;
    else
        val = opts.(valname);
    end
end