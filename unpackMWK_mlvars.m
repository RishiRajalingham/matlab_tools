function unpackMWK_mlvars(filename)

    %% Initialize
    addpath('/Library/Application Support/MWorks/Scripting/Matlab')
    addpath('/mindhive/dicarlolab/u/rishir/lib/matlab_tools/MWK_MATLAB/')
    if nargin < 1
        [fn, pathfn] = uigetfile('.mwk');
        filename = [pathfn, fn];
    end
    outfn = [filename(1:end-3), 'mat'];
    
    codec = getCodecs(filename);
    codec = codec.codec;
    tag_ind = find(strcmp(fieldnames(codec), 'tagname'));
    code_ind = find(strcmp(fieldnames(codec), 'code'));
    codec = struct2cell(codec);
    codec = codec([code_ind,tag_ind],:);
    
    %% Get codes for event tagnames
    ncodes = 15;
    codesOI = nan(ncodes,1); 
    
    % Scalar variables
    t = strcmp('sync', codec(2,:));
    tmp = codec(1,t);
    syncCode = double(tmp{1,1});
        
    for var_id = 1:ncodes
        try
            varname = ['ml_v', num2str(var_id-1)];
            t = strcmp(varname, codec(2,:));
            tmp = codec(1,t);
            codesOI(var_id) = double(tmp{1,1});
        catch rr_excep
            display(rr_excep)
        end
    end
    
    ncodes = sum(~isnan(codesOI));
    display('Retrieved codes...');

    %% Get trial times to sync from
    [tsync, sync] = getTimeStampedData(syncCode);
    
    tcode = cell(ncodes,1);
    xcode = cell(ncodes,1);
    for ci = 1:ncodes
        [tcode{ci}, xcode{ci}] = getTimeStampedData(codesOI(ci));
    end
    
    t1 = find(sync == 1);   t0 = find(sync == 0);
    t10 = intersect(t1, t0-1); % all 1s followed by 0
    N = length(t10);
    display(['Collecting ', num2str(N), ' trials of scalar data ']);
    
    %% Get scalar event data from trials
    allTrials = nan(N,ncodes);
    for i = 1:N
        start_sync = double(tsync(t10(i)));
        end_sync = double(tsync(t10(i)+1));
        for ci = 1:ncodes
            t_ = start_sync < tcode{ci} & tcode{ci} < end_sync;
            allTrials(i,ci) = xcode{ci}(t_);
        end
    end
    save(outfn, 'allTrials');
    display('Saved trials');
    
    %% Helper functions

    function [t,x] = getTimeStampedData(eventcode)
        events = getEvents(filename, eventcode);
        data_i = find(strcmp(fieldnames(events), 'data'));
        time_i = strcmp(fieldnames(events), 'time_us');
        
        events_vec = struct2cell(events);
        t = double(cell2mat(events_vec(time_i,:)));
        try
            x = double(cell2mat(events_vec(data_i,:)));
        catch % in case of different data types on different trials (??)
            tmp = sq(events_vec(data_i,:,:));
            x = nan(length(tmp));
            for tmpi = 1:length(tmp)
                x(tmpi) = double(tmp{tmpi});
            end
        end
        
        
    end
    
end