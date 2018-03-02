function [choice_x, B, B_ci] = rr_fit_psychometric(signal, choice, numparams, signal_x)
% function [choice_x, B, B_ci] = rr_fit_psychometric(signal, choice, numparams, signal_x)

    t = ~isnan(signal(:)) & ~isnan(choice(:));
    signal = signal(t);
    choice = choice(t);
    LAMBDA_DEFAULT = 0;%1/6;
    GAMMA_DEFAULT = 1 - LAMBDA_DEFAULT;
    
    if ~exist('numparams', 'var'); numparams = 2; end;
    if ~exist('signal_x', 'var'); signal_x = signal; end;
    
    alpha_thres = 0.99;
    [choice_x, B, B_ci] = fit_psych(signal, choice, signal_x);
     
    function [choice_x, B, B_ci] = fit_psych(signal, choice, signal_x)
        bl = [-Inf -Inf, 0, 0];
        bu = [Inf Inf, 1, 1];
        opts = optimset('Display','off', 'MaxFunEvals', 1000);
        beta0 = [0 1,1, 0];
        
        [B,resnorm, residual, exitflag, output, lambda, jacobian] = ...
            lsqcurvefit(@psychometric_fn, beta0(1:numparams), signal, choice, ...
            bl(1:numparams), bu(1:numparams), opts);
        choice_x = psychometric_fn(B,signal_x);
        B_ci = nlparci(B, residual, 'jacobian', jacobian);
%         
%         [B,R,~,covb] = nlinfit(signal, choice, @psychometric_fn, beta0(1:numparams));
        
%         B_ci = [];
        
       
%         B_ci = nlparci(B,R,'Covar', covb, 'alpha', alpha_thres);
    end

    function [choice_x, B, B_ci] = fit_psych_v2(signal, choice, signal_x)
        beta0 = [0 1,1, 0];
        [B,R,~,covb] = nlinfit(signal, choice, @psychometric_fn_v2, beta0(1:numparams));
        choice_x = psychometric_fn_v2(B,signal_x);
        B_ci = nlparci(B,R,'Covar', covb, 'alpha', alpha_thres);
    end

    warning('off','all');

    function z = psychometric_fn(beta0,t_)
        alpha_ = beta0(1);
        beta_ = beta0(2);
        if numparams > 2; gamma_ = beta0(3); else gamma_ = GAMMA_DEFAULT; end;
        if numparams > 3; lambda_ = beta0(4); else lambda_ = 1-gamma_; end;
        z = lambda_ + gamma_ .* (1 + exp(-(alpha_ + beta_*t_))).^-1;
    end

    function z = psychometric_fn_v2(beta0,t_)
        alpha_ = beta0(1);
        beta_ = beta0(2);
        if numparams > 2; gamma_ = beta0(3); else gamma_ = GAMMA_DEFAULT; end;
        if numparams > 3; lambda_ = beta0(4); else lambda_ = 1-gamma_; end;
        exparg = -beta_ .* (t_ - alpha_);
        
        z = lambda_ + gamma_ .* (1 + exp(exparg).^-1);
    end

end