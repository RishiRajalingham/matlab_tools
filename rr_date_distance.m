function dd = rr_date_distance(d1, d2)
% format dates as yyyymmdd, distance in actual days between

    d1 = date_to_string(d1);
    d2 = date_to_string(d2);
    dd = daysact(d1, d2);

    function D = date_to_string(d)
        D = sprintf('%d/%d/%d', floor(d./10000), floor(mod(d,10000)./100), mod(d,100));
    end

end