function [choice_x, B, B_ci] = rr_fit_psychometrics_paired(signal, choice, indicator, numparams, signal_x)
% function [choice_x, B, B_ci] = rr_fit_psychometrics_paired(signal, choice, indicator, numparams, signal_x)
% given choice for 2 conditions (selected by indicator), find optimal set
% of parameters, where condition is modeled as a shift in the psych curve.

    warning('off','all');
    signal = signal(:);
    choice = choice(:);
    indicator = indicator(:);
    
    if ~exist('numparams', 'var') || isempty(numparams)
        numparams = 3; 
    end
    if ~exist('signal_x', 'var') || isempty(signal_x)
        signal_x = [signal, indicator]; 
    end
    
    [choice_x, B, B_ci] = fit_psych(signal, choice, indicator, signal_x);
     
    function [choice_x, B, B_ci] = fit_psych(signal, choice, indicator, signal_x)
        bl = [-Inf -Inf, -Inf, 0, 0];
        bu = [Inf Inf, Inf, 1, 1];
        opts = optimset('Display','off', 'MaxFunEvals', 1000);
        beta0 = [0 1,1,1, 0];
        
        B = lsqcurvefit(@psychometric_fn, beta0(1:numparams), ...
            [signal, indicator], choice, bl(1:numparams), bu(1:numparams), opts);
        choice_x = psychometric_fn(B, signal_x);
        B_ci = [];
    end


    function z = psychometric_fn(beta0,xdata)
        alpha_ = beta0(1);
        beta_ = beta0(2);
        shift_ = beta0(3);
        
        t_ = xdata(:,1);
        i_ = xdata(:,2);
        
        if numparams > 3; gamma_ = beta0(4); else gamma_ = 1; end;
        if numparams > 4; lambda_ = beta0(5); else lambda_ = 0; end;
        exparg = alpha_ + beta_.*t_ + shift_*i_;
        z = lambda_ + gamma_ .* (1 + exp(-exparg) ).^-1;
    end

  

end