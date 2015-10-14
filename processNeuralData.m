function processNeuralData(fn)

    if iscell(fn)
        for i = 1:length(fn)
            process_file(fn);
        end
    else
        process_file(fn);
    end
        
        
    function process_file(fn_)
        load(fn_, 'allTrials');
        if size(allTrials,2) < 600; return; end;
        rasters = allTrials(:,2:end);
        imid = allTrials(:,1);
        t = ~isnan(imid); 
        imid = imid(t); rasters = rasters(t,:);

        if ~isempty(strfind(fn, 'OBJ'))
            nImgPerObj = 100;
        elseif ~isempty(strfind(fn, 'FBOP'))
            nImgPerObj = 20;
        end
        objid = floor(imid ./ nImgPerObj) + 1;
        time_ax = -100-25+1:500-25;

%         if ~isempty(strfind(fn_, 'OBJv2')); objid = objid + 25; end;

        winconv = normpdf(-20:20, 0, 10);
        winconv = winconv ./ sum(winconv);

        binconv = ones(1,100);
        psth = conv2(rasters, winconv, 'same');
        FR = sum(rasters(:,200:350), 2);
        FRb = conv2(rasters, binconv, 'same');
        FRb = downsample(FRb', 50)';

        save(fn_, 'psth', 'time_ax', 'FR', 'FRb', 'objid', '-append');
        display(['Saved to ',  fn_]);
    end
end
