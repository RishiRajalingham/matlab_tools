function unpack_PlxMwk(filename_mwk, filename_plx, protocol)
% function unpack_PlxMwk(filename_mwk, filename_plx)
% rishir 

    imgset = 1; % obj
    if ~exist('protocol', 'var')
        protocol = 1; % 1-rsvp, 0-m2s
    end
    
    %% Init
    if nargin < 1
        [fn, path] = uigetfile('.mat', 'Select MWK file.');
        filename_mwk = [path, fn];
        
        [fn, path] = uigetfile('.mat', 'Select PLX file.');
        filename_plx = [path, fn];
    end
    
    
    %% Unpack data
    display('***unpack_PlxMwk***')
%     outfn = strrep(filename_plx, '.mat', '_FULL.mat');
    tmpname = strrep(filename_mwk, '_', '');
    outfn = strrep(tmpname, '.mat', '_FULL.mat');
    load(filename_mwk, 'strobes', 'trialSummary',...
        'rsvpFBOP_trials', 'rsvpOBJ_trials', 'm2s_trials', ...
        'blockFBOP', 'blockOBJ');
    load(filename_plx, 'M2S', 'RSVP');
    
    if protocol == 0
        trial_summary = trialSummary(m2s_trials,:);
        mwk_strobes = strobes(m2s_trials,:);
        plx_strobes = M2S.strobes;
        spikedata = M2S.spikedata;
    elseif imgset == 0; 
        trial_summary = trialSummary(rsvpFBOP_trials,:);
        mwk_strobes = strobes(rsvpFBOP_trials,:);
        plx_strobes = RSVP.strobes;
        spikedata = RSVP.spikedata;
    else
        trial_summary = trialSummary(rsvpOBJ_trials,:);
        mwk_strobes = strobes(rsvpOBJ_trials,:);
        plx_strobes = RSVP.strobes;
        spikedata = RSVP.spikedata;
    end
        
    %% Align data
    % N.B: plex and mwk data should contain the same number of trials but
    % sometimes closing the mwk stream gets rid of the last trial. 
    
    [pk, mk] = align_plx_mwk(plx_strobes, mwk_strobes);
    trial_summary = trial_summary(mk,:);
    units = fieldnames(spikedata);
    for ui = 1:length(units)
        spikedata.(units{ui}).FR = spikedata.(units{ui}).FR(pk,:,:);
    end
    trial_id = [plx_strobes(pk,2), mwk_strobes(mk,1)];
    
%     %% Sanity-checks
%     assert(length(M2S.trialID) == size(M2S_trialSummary,1));
%     assert(length(RSVP.trialID) == size(RSVP_trialSummary,1));
%     
%     assert(norm(M2S.trialID - M2S_trialSummary(:,end)) == 0);
%     assert(norm(RSVP.trialID - RSVP_trialSummary(:,end)) == 0);
%     
    %% Align RSVP responses for population
    if protocol == 1
        display('Aligning plx/mwk for RSVP')
        trial_summary = trial_summary(:,1:5);    
        units = fieldnames(spikedata);
        nUnits = length(units);
        [nRSVPTrials, nRSVPIms] = size(trial_summary);
        
        allTrials = cell(nUnits,1);
        for ui = 1:nUnits
            un = units{ui};
            allTrials{ui} = [];
            if str2double(un(strfind(un,'u')+1:end)) == 0;
                continue; 
            end;
        
            FR = spikedata.(units{ui}).FR;
            [aa,bb,nTimeBins] = size(FR);
            assert(aa == nRSVPTrials);
            assert(bb == nRSVPIms);

            im = reshape(trial_summary, nRSVPTrials*nRSVPIms, 1);
            n = reshape(FR, nRSVPTrials*nRSVPIms, nTimeBins);
            allTrials{ui} = [im, n];

%             [mu,im_i] = grpstats(n, im, {'mean', 'gname'});
%             im_i = str2double(im_i);
%             imResp(im_i, ui,:) = mu; 
        end

        tmp.trial_summary = trial_summary;
        tmp.allTrials = allTrials;
        tmp.trial_id = trial_id;
        tmp.units = units;
        if imgset == 0;
            RSVP_FBOP = tmp;
            if exist(outfn, 'file')
                save(outfn, 'RSVP_FBOP', 'filename_mwk', 'filename_plx', '-append');
            else
                save(outfn, 'RSVP_FBOP', 'filename_mwk', 'filename_plx');
            end
        else
            RSVP_OBJ = tmp;
            if exist(outfn, 'file')
                save(outfn, 'RSVP_OBJ', 'filename_mwk', 'filename_plx', '-append');
            else
                save(outfn, 'RSVP_OBJ', 'filename_mwk', 'filename_plx');
            end
        end
    end
    
    %% Align M2S responses for population
    if protocol == 0
        display('Aligning plx/mwk for M2S')
        trial_summary = trial_summary(:,1:6);    
        units = fieldnames(spikedata);
        nUnits = length(units);
        nM2STrials = size(trial_summary,1);
        
        allTrials = cell(nUnits,1);
        for ui = 1:nUnits
            un = units{ui};
            allTrials{ui} = [];
            if str2double(un(strfind(un,'u')+1:end)) == 0;
                continue; 
            end;
        
            FR = sq(spikedata.(units{ui}).FR);
            assert(size(FR,1) == nM2STrials);

            im = trial_summary(:,1);
            allTrials{ui} = [im, FR];
        end

        M2S.trial_summary = trial_summary;
        M2S.allTrials = allTrials;
        M2S.trial_id = trial_id;
        M2S.units = units;
        beep;
        
        if exist(outfn, 'file')
            save(outfn, 'M2S', 'filename_mwk', 'filename_plx', '-append');
        else
            save(outfn, 'M2S', 'filename_mwk', 'filename_plx');
        end

%         
%         M2S.trialSummary = M2S_trialSummary(:,1:4);
%         imResp = nan(2500, nUnits, 500);
% 
%         for ui = 1:nUnits
%             FR = sq(M2S.spikedata.(units{ui}).FR);
%             im = M2S.trialSummary(:,1);
% 
%             [mu,im_i] = grpstats(FR, im, {'mean', 'gname'});
%             im_i = str2double(im_i);
%             imResp(im_i, ui,:) = mu; 
%         end
% 
%         M2S.imResp = imResp;
    end
    
    display('Done saving.')
    %% Helper functions
    function [pk, mk] = align_plx_mwk(plx_strobes, mwk_strobes)
        SAMPLEOFF = 4;
        
        [plx_trial_keep, ~] = find(plx_strobes == SAMPLEOFF);
        plx_trial_keep = unique(plx_trial_keep);
        plx_trial_id = plx_strobes(plx_trial_keep,2);
        
        [mwk_trial_keep,~] = find(mwk_strobes == SAMPLEOFF);
        mwk_trial_keep = unique(mwk_trial_keep);
        mwk_trial_id = mwk_strobes(mwk_trial_keep,1);

        L_m = length(mwk_trial_id);
        L_p = length(plx_trial_id);
%         assert(L_m <= L_p); % plx recording contains all trials

        if L_m <= L_p % plx recording contains all trials
            [a,b] = xcorr(plx_trial_id, mwk_trial_id);
            [~,mi] = max(a);
            trial_subset = b(mi)+1 : b(mi) + L_m;
            plx_trial_keep = plx_trial_keep(trial_subset);

            plx_trial_id = plx_trial_id(trial_subset);
            assert(norm(plx_trial_id-mwk_trial_id) == 0);
        else
            [a,b] = xcorr(mwk_trial_id, plx_trial_id);
            [~,mi] = max(a);
            trial_subset = b(mi)+1 : b(mi) + L_p;
            
            if min(trial_subset) < 0;%plx contains trial id that are not in mwk (fbop).
                t = trial_subset > 0;
                plx_trial_keep = plx_trial_keep(t);
                plx_trial_id = plx_trial_id(t);
                trial_subset = trial_subset(t);
            end
            
            
            mwk_trial_keep = mwk_trial_keep(trial_subset);
            mwk_trial_id = mwk_trial_id(trial_subset);
            
            t = plx_trial_id == mwk_trial_id;
            plx_trial_id = plx_trial_id(t);
            mwk_trial_id = mwk_trial_id(t);
            
            assert(norm(plx_trial_id-mwk_trial_id) == 0);
        end
        
        pk = plx_trial_keep;
        mk = mwk_trial_keep;

    end

%% to run-after
% d = dir('plex_data/*.mat');
% for di = 1:length(d)
%     fn = ['plex_data/', d(di).name];
%     unpack_PlxMwk('MantoRSVP_20150110.mat', fn)
% end
    
end

