function unpackPLXData(inputarg)
% function unpackPLXData(inputarg)
% rishir - 2014/07
% Unpack spike rasters aligned to stimulus onset for match-to-sample and
% rsvp protocols. Also save trial ID for sync-ing with behavioural data and
% strobes/strobeTimes for validation.
% input is filename or cell of filenames

    %% init
    global filename 
    global TRIALSTART FIXACQUIRE SAMPLEON SAMPLEOFF TESTON SELECTIONMADE NOSELECTION RSVPEND FIXBROKEN
    global strobes strobe_times tscounts trial_strobes trial_strobeTimes
    
    %% Strobe library
      TRIALSTART = 1;
    FIXACQUIRE = 2;
    SAMPLEON = 3;
    SAMPLEOFF = 4;
    TESTON = 5;
    SELECTIONMADE = 6;
    NOSELECTION = 7;
    RSVPEND = 8;
%     FIXBROKEN = 9;
    FIXBROKEN = 0;
      
    %% Parse and run 
    display('***UnpackPlx***')
    addpath(genpath('/Volumes/dicarlo/u/rishir/lib/plexlib/'));
    if nargin < 1
        [fn, path] = uigetfile('.plx');
        filename = [path, fn];
    elseif iscell(inputarg)
        for ai = 1:length(inputarg)
            filename = inputarg{ai};
            extract_and_save();
            clear filename;
        end
    else
        filename = inputarg;
        extract_and_save();
        clear filename;
    end
    
  
    %% Helper functions
    
    function extract_and_save()
        % Extract info
        [strobes, strobe_times, tscounts] = get_info();
        run_sanity_check();

        [trial_strobes, trial_strobeTimes] = getM2S_trialstrobes();
        M2S = struct;
        M2S.trialID = trial_strobes(:,3);
        M2S.strobes = trial_strobes;
        M2S.strobe_times = trial_strobeTimes;
        M2S.spikedata = getSpikeData('M2S');

        [trial_strobes, trial_strobeTimes] = getRSVP_trialstrobes();
        RSVP = struct;
        RSVP.trialID = trial_strobes(:,2);
        RSVP.strobes = trial_strobes;
        RSVP.strobe_times = trial_strobeTimes;
        RSVP.spikedata = getSpikeData('RSVP');

        % Save data
        outfn = strrep(filename, '.plx', '.mat');
        save(outfn, 'M2S', 'RSVP', 'strobes', 'strobe_times');
        display(['Saved to ...', filename]);
    end

    function [strobes, strobe_times, tscounts] = get_info()
        % Returns all strobes and strobetimes, after correcting alignment,
        % converting to ms. Also return spike
        % count per channel/unit.
        [~,strobe_times,strobes] = plx_event_ts(filename, 257);
        tscounts = plx_info(filename, 0);
            % correct for plexon offset bug, convert to ms
        strobes = strobes(2:end); 
        strobe_times = strobe_times(1:end-1) .* 1000;
        
            % repeat strobes bug?? if the same strobe is sent twice -
            % ignore second
        t = find((strobes(1:end-1) == strobes(2:end)) & (strobes(1:end-1)~=0));
        t = setdiff(1:length(strobes), t);
        strobes = strobes(t);
        strobe_times = strobe_times(t);
        
    end

    function run_sanity_check()
        % check stimulus presentation is 100ms
        
        [~, sample_pres_t] = getStrobeSequence(SAMPLEON,SAMPLEOFF,1);
        assert(abs(mean(sample_pres_t) - 100) < 20);
    end

    function [s12, del_t] = getStrobeSequence(strobe1, strobe2, interval)
        % Return all instances of a strobe sequence (s1,s2) separated by
        % defined interval. Also return the time difference in ms.
        
        s1 = find(strobes == strobe1);
        t = (s1+interval > 0) & (s1+interval < length(strobes));
        s1 = s1(t);
        s12 = s1(strobes(s1+interval) == strobe2);
        del_t = strobe_times(s12+interval) - strobe_times(s12); 
    end

    %% Extraction functions
    
    function spikedata = getSpikeData(protocol)
        % Return all spike data for all channels with spikes. Structure
        % contains rasters and mean waveform. For M2S.
        spikedata = struct();
        for unit_i = 1:size(tscounts,1)
            for chan_i = 1:size(tscounts,2)
                if tscounts(unit_i, chan_i) == 0; continue; end;
                ch = chan_i - 1; u = unit_i - 1;
                if strcmp(protocol, 'M2S')
                    win = 600;
                    numwin = 1;
                elseif strcmp(protocol, 'RSVP')
                    win = 600;
                    numwin = 5;
                end
                [spks, meanwave] = get_data_unit(ch, u, win, numwin);
                spikedata.(['ch',num2str(ch), 'u', num2str(u)]).FR = spks;
                spikedata.(['ch',num2str(ch), 'u', num2str(u)]).waves = meanwave;
%                 spikedata.(['ch',num2str(ch), 'u', num2str(u)]).spkt = spkt;
            end
        end
    end

    function [trialSpikes, meanWave, allSpikes] = get_data_unit(ch, u, win, numwin)
        % Returns a raster of fixed window-size of spike counts for a
        % particular neuron (specified by channel/unit). For M2S.
        nTrials = size(trial_strobes,1); 
        trialSpikes = nan(nTrials, numwin, win);
        
        [~, spk_t] = plx_ts(filename, ch, u);
        spk_t = spk_t .* 1000;
%         allSpikes = spk_t;
        
        for t_i = 1:nTrials
            t_son = find(trial_strobes(t_i,:) == SAMPLEON);
            for s_i = 1:length(t_son)
                t_start = trial_strobeTimes(t_i,t_son(s_i)) -100;
                t_win = spk_t >= t_start & spk_t < t_start+win;
                spk_t_curr = spk_t(t_win) - t_start;
                trialSpikes(t_i,s_i,:) = hist(spk_t_curr, 0:win-1);
            end
        end 
        
        [~,~,~,w] = plx_waves_v(filename, ch, u);
        meanWave = mean(w);
    end

    function [trial_strobes, trial_strobeTimes] = getM2S_trialstrobes()
        % Return strobes of interest per trial for all completed trials,
        % and their corresponding strobe times. For M2S data.
        nStrobesPerTrial = 7;
        trial_end = getStrobeSequence(SELECTIONMADE,TESTON,-1);
        trial_end = trial_end(trial_end >= nStrobesPerTrial);
    
        nTrials = length(trial_end);
        fprintf(1, 'Extracting %d trials of M2S ... \n', nTrials);
        trial_strobes = nan(nTrials,nStrobesPerTrial);
        trial_strobeTimes = nan(nTrials,nStrobesPerTrial);

        for i = 1:nTrials
            win = (trial_end(i)-nStrobesPerTrial+1):trial_end(i);
            trial_strobes(i,:) = strobes(win);
            trial_strobeTimes(i,:) = strobe_times(win);
        end
        % sanity check
        if nTrials > 0
            ind = setdiff(1:nStrobesPerTrial, 3);
            assert(sum(nanstd(trial_strobes(:,ind))) == 0);
        end
    end

    function [trial_strobes, trial_strobeTimes] = getRSVP_trialstrobes()
        % Return strobes of interest per trial for all completed trials,
        % and their corresponding strobe times. For RSVP data.
        nStrobesPerTrial = 13;
        trial_start_i = find(strobes == TRIALSTART);
        trial_end_1 = getStrobeSequence(RSVPEND,SAMPLEOFF,-1);
        trial_end_2 = getStrobeSequence(FIXBROKEN,SAMPLEOFF,-1);
        trial_end_3 = getStrobeSequence(FIXBROKEN,SAMPLEON,-1);
        trial_end = [trial_end_1; trial_end_2; trial_end_3];
        trial_end = trial_end(trial_end >= nStrobesPerTrial);
        trial_end = sort(trial_end);
    
        nTrials = length(trial_end);
        fprintf(1, 'Extracting %d trials of RSVP ... \n', nTrials);
        trial_strobes = nan(nTrials,nStrobesPerTrial);
        trial_strobeTimes = nan(nTrials,nStrobesPerTrial);
%         keep_trial = false(nTrials,1);
        for i = 1:nTrials
            tmp_start = find(trial_start_i < trial_end(i),1, 'last');
            trial_start = trial_start_i(tmp_start);
            win = trial_start:trial_end(i);
%             keep_trial(i) = ismember(SAMPLEOFF, strobes(win));
            trial_strobes(i,1:length(win)) = strobes(win);
            trial_strobeTimes(i,1:length(win)) = strobe_times(win);
        end
%         
        keep_trial = true(nTrials,1);
        
        % sanity check
%         if nTrials > 0
%             ind = setdiff(1:nStrobesPerTrial, 2);
%             assert(sum(nanstd(trial_strobes(:,ind))) == 0);
%         end
        trial_strobes = trial_strobes(keep_trial,:);
        trial_strobeTimes = trial_strobeTimes(keep_trial,:);
        fprintf(1, 'Keeping %d trials of RSVP ... \n', size(trial_strobes,1));
    end

    
end