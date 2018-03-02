 function D = rr_smooth2d(inputX, sigma, winpxl, nanopt)
    % 2d smoothing with gaussian of width sigma
    if ~exist('nanopt', 'var')
        nanopt = 0; %default nan options
    end
    if ~exist('winpxl', 'var')
        winpxl = 1;
    end
    if isnan(sigma)
        D = inputX;
    else
        [m,n] = size(inputX);
        onesX = ones(size(inputX)); 
        onesX(isnan(inputX)) = nan;
        
        mu = [0 0];
        Sigma = [sigma.^2 0; 0 sigma.^2];
        x1 = -winpxl:winpxl; 
        x2 = -winpxl:winpxl; 
        
        [X1,X2] = meshgrid(x1,x2);
        F = mvnpdf([X1(:) X2(:)],mu,Sigma);
        F = reshape(F,length(x2),length(x1));

        if nanopt == 0
            D = nanconv(inputX, F);
            D2 = nanconv(onesX, F);
        elseif nanopt == 1
            D = nanconv(inputX, F, 'noedge', 'nanout');
            D2 = nanconv(onesX, F, 'noedge', 'nanout');
        end
        D = D ./ D2;
    end
    
end