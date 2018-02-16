function unpackMWK_mlvars(filename, codenames, codenames_analog)
% unpackMWK_mlvars(filename, codenames, codenames_analog)


%% Opt
get_dat = true;
get_analog = true;
get_spikes = true;
get_hashes = false;
check_led_error = false;
allTrials = []; 
allTrials.codenames = codenames;
allTrials.codenames_analog = codenames_analog;

%% Initialize
addpath('/Library/Application Support/MWorks/Scripting/Matlab')
addpath(genpath('/mindhive/dicarlolab/u/rishir/lib/matlab_tools/MWK_MATLAB/'))
if nargin < 1
    [fn, pathfn] = uigetfile('.mwk');
    filename = [pathfn, fn];
end
outfn = strrep(filename, '.mwk', '.mat');

%% Get codenames and codes
codec = getCodecs(filename);
codec = codec.codec;
tag_ind = find(strcmp(fieldnames(codec), 'tagname'));
code_ind = find(strcmp(fieldnames(codec), 'code'));
codec = struct2cell(codec);
codec = codec([code_ind,tag_ind],:);
allcodecs = codec(2,:);

if ~exist('codenames', 'var')
    codenames = get_default_codenames();
end
ncodes = length(codenames);

if ~exist('codenames_analog', 'var')
    codenames_analog = {};
end
ncodes_analog = length(codenames_analog);

%% Get codes for event tagnames
codesOI = nan(ncodes,1);

% Scalar variables
t = strcmp('sync', codec(2,:));
tmp = codec(1,t);
syncCode = double(tmp{1,1});

t = strcmp('#announceMessage', codec(2,:));
tmp = codec(1,t);
messageCode = double(tmp{1,1});

t = strcmp('#announceStimulus', codec(2,:));
tmp = codec(1,t);
stimulisCode = double(tmp{1,1});

t = strcmp('#stimDisplayUpdate', codec(2,:));
tmp = codec(1,t);
displayCode = double(tmp{1,1});

t = strcmp('plex_spike', codec(2,:));
tmp = codec(1,t);
spikeCode = double(tmp{1,1});

for var_id = 1:ncodes
    try
        varname = codenames{var_id};
        t = strcmp(varname, codec(2,:));
        tmp = codec(1,t);
        codesOI(var_id) = double(tmp{1,1});
    catch
        fprintf(1, 'Exception for code %s \n',varname);
    end
end
fprintf(1, '%s : Retrieved %d codes... \n', filename, ncodes);

analog_codesOI = nan(ncodes_analog,1);
for var_id = 1:ncodes_analog
    try
        varname = codenames_analog{var_id};
        t = strcmp(varname, codec(2,:));
        tmp = codec(1,t);
        analog_codesOI(var_id) = double(tmp{1,1});
    catch
        fprintf(1, 'Exception for code \n %s',varname);
    end
end

%% Get trial times to sync from
[tsync, sync] = getTimeStampedData(syncCode);

tcode = cell(ncodes,1);
xcode = cell(ncodes,1);
for ci = 1:ncodes
    if ~isnan(codesOI(ci))
        [tcode{ci}, xcode{ci}] = getTimeStampedData(codesOI(ci));
    end
end

tcode_analog = cell(ncodes_analog,1);
xcode_analog = cell(ncodes_analog,1);
for ci = 1:ncodes_analog
    if ~isnan(analog_codesOI(ci))
        [tcode_analog{ci}, xcode_analog{ci}] = getTimeStampedData(analog_codesOI(ci));
    end
end

[start_sync_t, end_sync_t] = get_sync_times(xcode, tcode, sync, tsync);

N = length(start_sync_t);
fprintf(1, 'Collecting %d trials of scalar data \n', N);

ttrial_mu = nanmean(end_sync_t - start_sync_t);
fprintf(1, 'Avg trial duration : %f seconds \n', ttrial_mu./10^3);
beep;
%% Get scalar event data from trials
if get_dat
    nvals = 50; % keep track of 5 first values of code-- be careful to use last
    allTrials.dat = nan(N,ncodes, nvals);
    allTrials.time = nan(N,ncodes, nvals);
    for i = 1:N
        start_sync = start_sync_t(i);
        end_sync = end_sync_t(i);
        for ci = 1:ncodes
            t_ = find(start_sync < tcode{ci} & tcode{ci} <= end_sync, nvals, 'first');
            if isempty(t_)
                continue;
            end
            allTrials.dat(i,ci, 1:length(t_)) = xcode{ci}(t_);
            allTrials.time(i,ci, 1:length(t_)) = tcode{ci}(t_);
        end
    end
end

%% Get analog variables from trials
if get_analog
    analog_duration = 2000;
    allTrials.analog = cell(ncodes_analog,1);
    if ncodes_analog > 0
        for ci = 1:ncodes_analog
            allTrials.analog{ci} = nan(N,analog_duration);
        end
        for i = 1:N
            start_sync = start_sync_t(i);
            end_sync = end_sync_t(i);
            for ci = 1:ncodes_analog
                t_ = find(start_sync < tcode_analog{ci} & tcode_analog{ci} < end_sync);
                tt = tcode_analog{ci}(t_); xx = xcode_analog{ci}(t_);

                if isempty(t_)
                    continue; 
                end
                try
                    t_align = start_sync:1:end_sync;
                    tlen = min(analog_duration, length(t_align));
                    t_align = t_align(1:tlen);
                    allTrials.analog{ci}(i,1:tlen) = interp1(tt,xx,t_align);
                catch
                    fprintf(1, 'exception analog align \n')
                end
            end
        end
    end
end

%% Get spike data from trials
if get_spikes
%     [t_disp, x_disp] = getTimeStampedData(displayCode);
    [t_spk, ~] = getTimeStampedData(spikeCode);    
    nspikevals = 10000;
    allTrials.spikes.time = nan(N,nspikevals);
    for i = 1:N
        start_sync = start_sync_t(i);
        end_sync = end_sync_t(i);    
        t_ = find(start_sync < t_spk & t_spk <= end_sync, nspikevals, 'first');
        allTrials.spikes.time(i,1:length(t_)) = t_spk(t_);
        
    end
end

%% get hashes data from trials
if get_hashes
    allTrials.hash = nan(N,40);
    [tim,xim] = getStimulusMessages(stimulisCode);
    for i = 1:N
        start_sync = start_sync_t(i);
        end_sync = end_sync_t(i);
        t_ = find(start_sync < tim & tim <= end_sync, 1, 'first');
        if isempty(t_)
            continue; 
        end
        allTrials.hash(i,:) = xim(t_,:);
    end
end

%% Get announced messages from trials

if check_led_error
    good_trials = nan(N,1);
    error_status = 4;
    [tmess,~] = getErrorMessages(messageCode, error_status);
    
    for i = 1:N
        start_sync = start_sync_t(i);
        end_sync = end_sync_t(i);
        t_ = sum(start_sync < tmess & tmess < end_sync);
        good_trials(i) = t_ == 0;
    end
    good_trials = logical(good_trials);
    allTrials.good_trials = good_trials;
end


%%

save(outfn, 'allTrials', 'codenames',  'codenames_analog');
fprintf(1, 'Saved to %s \n', outfn);

%% Helper functions

    function codenames = get_default_codenames()
        ncodes_ = 0;
        for i_ = 1:length(allcodecs)
            if ~isempty(strfind(allcodecs{i_}, 'ml_v'))
                ncodes_ = ncodes_ + 1;
            end
        end
        codenames = cell(ncodes_,1);
        for var_id_ = 1:ncodes_
            codenames{var_id_} = ['ml_v', num2str(var_id_-1)];
        end
    end

    function [t,x] = getTimeStampedData(eventcode)
        events = getEvents(filename, eventcode);
        data_i = find(strcmp(fieldnames(events), 'data'));
        time_i = strcmp(fieldnames(events), 'time_us');
        
        events_vec = struct2cell(events);
        t = double(cell2mat(events_vec(time_i,:))) ./ 1000; % in ms
        try
            x = double(cell2mat(events_vec(data_i,:)));
        catch % in case of different data types on different trials (??)
            fprintf(1, 'Exception for getTimeStampedData %d \n', eventcode);
            tmp_ = sq(events_vec(data_i,:,:));
            x = nan(length(tmp_),1);
            for tmpi = 1:length(tmp_)
                x(tmpi) = double(tmp_{tmpi});
            end
        end
    end

    function [terr,xerr] = getErrorMessages(eventcode, error_status)
        events = getEvents(filename, eventcode);
        xerr = []; terr = [];
        error_message = ['status: ', num2str(error_status)];
        for ei = 1:length(events)
            if ~isstruct(events(ei).data); continue; end;
            i_ = strfind(events(ei).data.message, 'ERROR');
            j_ = strfind(events(ei).data.message, error_message);
            if ~isempty(i_) && ~isempty(j_)
                xerr = cat(1, xerr, error_status);
                terr = cat(1, terr, events(ei).time_us);
            end
        end
    end

    function [tim,xim] = getStimulusMessages(eventcode)
        events = getEvents(filename, eventcode);
        xim = []; tim = [];
        for ei = 1:length(events)
            if ~isstruct(events(ei).data); continue; end;
            if ismember('file_hash', fieldnames(events(ei).data))
                xim = cat(1, xim, events(ei).data.file_hash);
                tim = cat(1, tim, events(ei).time_us);
            end
        end
    end

    function [start_sync_t, end_sync_t] = get_sync_times(xcode, tcode, sync, tsync)
        if true
            t1 = find(sync == 1);   t0 = find(sync == 0);
            t10 = intersect(t1, t0-1); % all 1s followed by 0
            start_sync_t = tsync(t10);
            end_sync_t = tsync(t10+1);
        else % local variable sync bug
            t_tmp = find(diff(xcode{15}) ~= 0); % repeat val of temp_pre
            nT = length(t_tmp);
            start_sync_t = nan(nT,1);
            end_sync_t = nan(nT,1);
            
            for ti = 1:length(t_tmp)
                if xcode{15}(t_tmp(ti)) == 0
                    continue;
                end
                start_sync_t(ti) = tcode{15}(t_tmp(ti));
                sync_end_i = find(tcode{16} > start_sync_t(ti), 1, 'first');
                end_sync_t(ti) = tcode{16}(sync_end_i);
            end
        end
        
    end


end