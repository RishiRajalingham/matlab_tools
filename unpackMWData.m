function unpackMWData(filename, getEyeData, getStrobeData, getSpikeData)
% function unpackMWData(filename, getEyeData, getStrobeData, getSpikeData)
% rishir 

    %% Initialize
    addpath('/Library/Application Support/MWorks/Scripting/Matlab')
    if nargin < 1
        [fn, pathfn] = uigetfile('.mwk');
%         [fn, pathfn] = uigetfile();
        filename = [pathfn, fn];
    end
    
%     outfn = [pathfn, strrep(fn, 'mwk', 'mat')];
    outfn = [filename(1:end-3), 'mat'];
    
    if ~exist('getEyeData', 'var') ; getEyeData = 0; end;
    if ~exist('getStrobeData', 'var'); getStrobeData = 1; end;
    if ~exist('getSpikeData', 'var'); getSpikeData = 1; end;
    if ~exist('getDiodeData', 'var'); getDiodeData = 1; end;
    
    %% Get codes for event tagnames
    codec = getCodecs(filename);
    codec = codec.codec;
    
    % CODE LIBRARY
    FIXBROKEN = 1;
    SUCCESS = 2;
    FAILURE = 3;
    IGNORE = 4;
    SAMPLEIMG = 5;
    TESTMATCHIMG = 6;
    TESTNONMATCHIMG = 7;
    IMGLOCATION = 8;
    SELECTIONOPEN = 9;
    EYEH = 10;
    EYEV = 11;
    EYEINSAMPLEWINDOW = 12;
    EYEINOTHERWINDOW = 13;
    STROBE = 14;
    PLEXSPK = 15;
    BLOCKFBOP = 16;
    BLOCKOBJ25 = 17;
    BLOCKOBJ25v2 = 19;
    VERTREFRESH= 18;
    
    syncCode = nan;
    codesOI = nan(17,1); 
      
    for i = 1:length(codec)
        if (strcmp(codec(i).tagname, 'sync'))
            syncCode = [codec(i).code];
        elseif (strcmp(codec(i).tagname, 'fixation_broken'))
            codesOI(FIXBROKEN) = [codec(i).code];
        elseif (strcmp(codec(i).tagname, 'success'))
            codesOI(SUCCESS) = [codec(i).code];
        elseif (strcmp(codec(i).tagname, 'failure'))
            codesOI(FAILURE) = [codec(i).code];
        elseif (strcmp(codec(i).tagname, 'ignore'))
            codesOI(IGNORE) = [codec(i).code];
        elseif (strcmp(codec(i).tagname, 'sample_image_index'))
            codesOI(SAMPLEIMG) = [codec(i).code];
        elseif (strcmp(codec(i).tagname, 'test_match_image_index'))
            codesOI(TESTMATCHIMG) = [codec(i).code];
        elseif (strcmp(codec(i).tagname, 'test_nonmatch_image_index'))
            codesOI(TESTNONMATCHIMG) = [codec(i).code];
        elseif (strcmp(codec(i).tagname, 'image_location_index'))
            codesOI(IMGLOCATION) = [codec(i).code];
        elseif (strcmp(codec(i).tagname, 'selection_open'))
            codesOI(SELECTIONOPEN) = [codec(i).code];
        elseif (strcmp(codec(i).tagname, 'eye_h'))
            codesOI(EYEH) = [codec(i).code];
        elseif (strcmp(codec(i).tagname, 'eye_v'))
            codesOI(EYEV) = [codec(i).code];
        elseif (strcmp(codec(i).tagname, 'eye_in_sample_image_window'))
            codesOI(EYEINSAMPLEWINDOW) = [codec(i).code];
        elseif (strcmp(codec(i).tagname, 'eye_in_other_image_window'))
            codesOI(EYEINOTHERWINDOW) = [codec(i).code];
        elseif (strcmp(codec(i).tagname, 'strobe'))
            codesOI(STROBE) = [codec(i).code];  
        elseif (strcmp(codec(i).tagname, 'plex_spike'))
            codesOI(PLEXSPK) = [codec(i).code];  
        elseif (strcmp(codec(i).tagname, 'num_blocks_shown_fbop'))
            codesOI(BLOCKFBOP) = [codec(i).code];  
        elseif (strcmp(codec(i).tagname, 'num_blocks_shown_obj25'))
            codesOI(BLOCKOBJ25) = [codec(i).code];  
        elseif (strcmp(codec(i).tagname, 'num_blocks_shown_obj25_v2'))
            codesOI(BLOCKOBJ25v2) = [codec(i).code];  
        elseif (strcmp(codec(i).tagname, 'vertical_refresh'))
            codesOI(VERTREFRESH) = [codec(i).code];  
        end       
    end
    
    display('Retrieved codes...');
    
    %% Get trial times to sync from
    syncEvents = getEvents(filename, syncCode);
    sync = zeros(length(syncEvents),1);
    for i = 1:length(syncEvents)
        sync(i) = syncEvents(i).data;
    end
    
    display(length(sync))
    display(sum(sync))
    
    t1 = find(sync == 1);   t0 = find(sync == 0);
    t10 = intersect(t1, t0-1); % all 1s followed by 0
    
    %% Declare variables to fill in with extraction
    N = length(t10);
    display(['Collecting ', num2str(N), ' trials of data '])
    
    trialSummary = nan(N,6);
    
    m2s_trials = false(N,1);
    rsvpFBOP_trials = false(N,1);
    rsvpOBJ_trials = false(N,1);
    
    blockFBOP = nan(N,1);
    blockOBJ = nan(N,1);
    
    strobes = nan(N,10);
    strobeTimes = nan(N,10);
    
    spikeTimes = nan(N,10000);
    
    diode = nan(N,10,2);
    
    analogEyeX = nan(N,5000);
    analogEyeY = nan(N,5000);
    
    triggerEyeSample = nan(N,5000);
    triggerEyeOther  = nan(N,5000);
    
    %% Get event data from trials
    
    for i = 1:N
        
        start_sync = double(syncEvents(t10(i)).time_us);
        end_sync = double(syncEvents(t10(i)+1).time_us);
        
        fbopEvent = getEvents(filename,codesOI(BLOCKFBOP), start_sync, end_sync);
        objEvent = getEvents(filename,codesOI(BLOCKOBJ25), start_sync, end_sync);
        try
            obj2Event = getEvents(filename,codesOI(BLOCKOBJ25v2), start_sync, end_sync);
        catch 
            obj2Event = [];
        end
        
        try
            selectionEvent = getEvents(filename,codesOI(SELECTIONOPEN), start_sync, end_sync);
        catch
            selectionEvent = [];
        end
        
        if ~isempty(fbopEvent); 
            rsvpFBOP_trials(i) = true; 
            [trialSummary(i,:), blockFBOP(i)] = unpackRSVP(start_sync, end_sync);
        elseif ~isempty(objEvent) 
            rsvpOBJ_trials(i) = true; 
            [trialSummary(i,:), blockOBJ(i)] = unpackRSVP(start_sync, end_sync);
        elseif ~isempty(obj2Event) 
            rsvpOBJ_trials(i) = true; 
            [trialSummary(i,:), blockOBJ(i)] = unpackRSVP(start_sync, end_sync);
        elseif ~isempty(selectionEvent); 
            m2s_trials(i) = true; 
            trialSummary(i,:) = unpackM2S(start_sync, end_sync);
        end;
        
%         if getEyeData == 1
%             analog_eyeh = getEvents(filename, codesOI(EYEH), start_sync, end_sync);
%             analog_eyev = getEvents(filename, codesOI(EYEV), start_sync, end_sync);
%             ttl_eyesample = getEvents(filename, codesOI(EYEINSAMPLEWINDOW), start_sync, end_sync);
%             ttl_eyeother = getEvents(filename, codesOI(EYEINOTHERWINDOW), start_sync, end_sync);
%         end
%         
        
        if getSpikeData == 1
            spikeEvents = getEvents(filename, codesOI(PLEXSPK), start_sync, end_sync);
            [spikeT,~] = getEventTimeSeries(spikeEvents);
            spikeTimes(i,1:length(spikeT)) = spikeT;
        end
        
        if getStrobeData == 1
            strobeEvents = getEvents(filename, codesOI(STROBE), start_sync, end_sync);
            for se_i = 1:length(strobeEvents)
                strobes(i,se_i) = strobeEvents(se_i).data;
                strobeTimes(i,se_i) = double(strobeEvents(se_i).time_us)/1000;
            end
        end
        
        if getDiodeData == 1
            diodeEvents = getEvents(filename, codesOI(VERTREFRESH), start_sync, end_sync);
            for se_i = 1:length(diodeEvents)
                diode(i,se_i,1) = diodeEvents(se_i).data;
                diode(i,se_i,2) = double(diodeEvents(se_i).time_us)/1000;
            end
        end
        
        if mod(i, floor(N/100)) == 1
            display(['Done trial ', num2str(i), ' of ', num2str(N)])
        end
    end
    
    M2S.trialSummary = [trialSummary(m2s_trials,:), strobes(m2s_trials,1)];
    RSVP_FBOP.trialSummary = [trialSummary(rsvpFBOP_trials,:), strobes(rsvpFBOP_trials,1)];
    RSVP_OBJ.trialSummary = [trialSummary(rsvpOBJ_trials,:), strobes(rsvpOBJ_trials,1)];
        
    
    save(outfn, 'trialSummary', 'm2s_trials', 'rsvpFBOP_trials', 'rsvpOBJ_trials', ...
        'analogEyeX', 'analogEyeY','triggerEyeSample','triggerEyeOther', ...
        'M2S', 'RSVP_FBOP', 'RSVP_OBJ', 'blockOBJ', 'blockFBOP', ...
        'strobes', 'strobeTimes', 'spikeTimes', 'diode');
    
    
    
    function trial_data = unpackM2S(start_sync, end_sync)
    % Extract summary of one trial of m2s data.
        trial_data = nan(1,6);
        successEvent = getEvents(filename,codesOI(SUCCESS), start_sync, end_sync);
        punishEvent = getEvents(filename,codesOI(FAILURE), start_sync, end_sync);
        ignoreEvent = getEvents(filename,codesOI(IGNORE), start_sync, end_sync);
        selecEvent = getEvents(filename,codesOI(SELECTIONOPEN), start_sync, end_sync);
               
        select_start_sync = selecEvent(1).time_us;
%         select_end_sync = selecEvent(2).time_us;
        
        if ~isempty(successEvent)
            trial_data(5) = (double(successEvent(1).time_us) - select_start_sync)/1000; % reaction time
            trial_data(6) = 1; % correct match
        elseif ~isempty(punishEvent)
            trial_data(5) = (double(punishEvent(1).time_us) - select_start_sync)/1000; % reaction time
            trial_data(6) = -1; % incorrect match
        elseif ~isempty(ignoreEvent)
            trial_data(5) = NaN; % reaction time
            trial_data(6) = 0; % no selection
        end  

        imgEvents_sampleImgIndex = getEvents(filename,codesOI(SAMPLEIMG), start_sync, end_sync);
            trial_data(1) = imgEvents_sampleImgIndex(1).data; % current sample index
        imgEvents_matchIndex = getEvents(filename,codesOI(TESTMATCHIMG), start_sync, end_sync);
            trial_data(2) = imgEvents_matchIndex(1).data; % current match index
        imgEvents_nonmatchIndex = getEvents(filename,codesOI(TESTNONMATCHIMG), start_sync, end_sync);
            trial_data(3) = imgEvents_nonmatchIndex(1).data; % current nonmatch index
        imgLocEvents_matchLocation = getEvents(filename,codesOI(IMGLOCATION), start_sync, end_sync);
            trial_data(4) = imgLocEvents_matchLocation(1).data; % match location
    end
    

    function [trial_data, block] = unpackRSVP(start_sync, end_sync)
    % Extract summary of one trial of rsvp data.
        trial_data = nan(1,6);        
        img_events = getEvents(filename,codesOI(SAMPLEIMG), start_sync, end_sync);
        fbop_event = getEvents(filename,codesOI(BLOCKFBOP), start_sync, end_sync);
        obj_event = getEvents(filename,codesOI(BLOCKOBJ25), start_sync, end_sync);
        try
            obj2_event = getEvents(filename,codesOI(BLOCKOBJ25v2), start_sync, end_sync);
        catch 
            obj2_event = [];
        end
        
        successEvent = getEvents(filename,codesOI(SUCCESS), start_sync, end_sync);
        fixEvent = getEvents(filename,codesOI(FIXBROKEN), start_sync, end_sync);
        
        for i_i = 1:length(img_events)
            trial_data(i_i) = img_events(i_i).data; % current sample index
        end
        
        if ~isempty(fbop_event); block = fbop_event(1).data; end; % block of fbop
        if ~isempty(obj_event); block = obj_event(1).data; end; % block of obj25  
        if ~isempty(obj2_event); block = obj2_event(1).data; end; % block of obj25  
        if ~isempty(successEvent); trial_data(6) = 1; end; % completed trial
        if ~isempty(fixEvent); trial_data(6) = -2; end; % broken fixation
        
    end
    
%     function getTTLData(codeoi, start_sync, end_sync)
%     end
%     
%     function analog_dat = getAnalogData(code_, start_sync, end_sync)
% 
%             analog_event = getEvents(filename,codesOI(code_), select_start_sync, select_end_sync);
%             [eyeh_time, eyeh_val] = getEventTimeSeries(analog_event);
% 
%             eyeSampleEvents = getEvents(filename,codesOI(EYEINSAMPLEWINDOW), select_start_sync, select_end_sync);
%             eyeOtherEvents = getEvents(filename,codesOI(EYEINOTHERWINDOW), select_start_sync, select_end_sync);
% 
%             select_start_sync_ms = double(select_start_sync)/1000;
%             select_end_sync_ms = double(select_end_sync)/1000;
%             time_resampled = select_start_sync_ms:1:select_end_sync_ms; %1ms bins
% 
%             [eyeh_time, eyeh_val] = getEventTimeSeries(eyehEvents);
%             [eyev_time, eyev_val] = getEventTimeSeries(eyevEvents);
% 
%             [eyeSample_time, eyeSample_val] = getEventTimeSeries(eyeSampleEvents);
%             [eyeOther_time, eyeOther_val] = getEventTimeSeries(eyeOtherEvents);
% 
%             tmp = interp1q(eyeh_time(:), eyeh_val(:), time_resampled(:))';
%                 analogEyeX(i,1:length(tmp)) = tmp;  % analog eye x position
%             tmp = interp1q(eyev_time(:), eyev_val(:), time_resampled(:))';
%                 analogEyeY(i,1:length(tmp)) = tmp; % analog eye y position
% 
%             triggerEyeSample(i,ceil(eyeSample_time - select_start_sync_ms)) = eyeSample_val;
%             triggerEyeOther(i,ceil(eyeOther_time - select_start_sync_ms)) = eyeOther_val;
%         end
% 
%     end
%     
    function [t,v] = getEventTimeSeries(ev)
    % Get data and time of single event
        nEv = length(ev);
        t = zeros(1,nEv);
        v = zeros(1,nEv);
        
        for evi = 1:length(ev)
            t(evi) = double(ev(evi).time_us)/1000;
            v(evi) = double(ev(evi).data);
        end
    end
    
end


%%
%                 selectionEvent = getEvents(filename,codesOI(SELECTIONOPEN), start_sync, end_sync);
%                 select_start_sync = selectionEvent(1).time_us;
%                 select_end_sync = selectionEvent(2).time_us;
% 
%                 if getEyeData == 1   % get analog eye and digital eye trigger 
%                     eyehEvents = getEvents(filename,codesOI(EYEH), select_start_sync, select_end_sync);
%                     eyevEvents = getEvents(filename,codesOI(EYEV), select_start_sync, select_end_sync);
% 
%                     eyeSampleEvents = getEvents(filename,codesOI(EYEINSAMPLEWINDOW), select_start_sync, select_end_sync);
%                     eyeOtherEvents = getEvents(filename,codesOI(EYEINOTHERWINDOW), select_start_sync, select_end_sync);
% 
%                     select_start_sync_ms = double(select_start_sync)/1000;
%                     select_end_sync_ms = double(select_end_sync)/1000;
%                     time_resampled = select_start_sync_ms:1:select_end_sync_ms; %1ms bins
% 
%                     [eyeh_time, eyeh_val] = getEventTimeSeries(eyehEvents);
%                     [eyev_time, eyev_val] = getEventTimeSeries(eyevEvents);
% 
%                     [eyeSample_time, eyeSample_val] = getEventTimeSeries(eyeSampleEvents);
%                     [eyeOther_time, eyeOther_val] = getEventTimeSeries(eyeOtherEvents);
% 
%                     tmp = interp1q(eyeh_time(:), eyeh_val(:), time_resampled(:))';
%                         analogEyeX(i,1:length(tmp)) = tmp;  % analog eye x position
%                     tmp = interp1q(eyev_time(:), eyev_val(:), time_resampled(:))';
%                         analogEyeY(i,1:length(tmp)) = tmp; % analog eye y position
% 
%                     triggerEyeSample(i,ceil(eyeSample_time - select_start_sync_ms)) = eyeSample_val;
%                     triggerEyeOther(i,ceil(eyeOther_time - select_start_sync_ms)) = eyeOther_val;
%                 end


                
            
            
        

