function TEST_rr_exact_test_decode()
% function TEST_rr_exact_test_decode()
% code to test the conceptual validity of statistical testing in high
% dimensions by first decoding.

    eff_size = [0.1,0.5,1,2,5];
    ndims = [1,5,10,100];
    niter = 10;

    p_vals_null = nan(niter, length(ndims), 2);
    p_vals_true = nan(niter, length(ndims), length(eff_size), 2);
    for iter = 1:niter
        for di = 1:length(ndims)
            d = ndims(di);
            X = randn(1000,d);
            X_2 = randn(1000,d);
            out_1 = rr_exact_test_decode(X,X_2);
            p_vals_null(iter,di,2) = min(rr_exact_test(X,X_2)).*d;
            p_vals_null(iter,di,1) = out_1.p;

            for fi = 1:length(eff_size)
                f = eff_size(fi);
                Y = randn(1000,d) + f;
                out_2 = rr_exact_test_decode(X,Y);
                p_vals_true(iter,di,di,2) = min(rr_exact_test(X,Y)).*d;
                p_vals_true(iter,di,fi) = out_2.p;
            end
        end

    end
    test_stats = [];
    test_stats.p_vals_null = p_vals_null;
    test_stats.p_vals_true = p_vals_true;
    
    save('test_stats_rr_exact_test_decode.mat', 'test_stats');

end