function processPlxUnits(datapath)
% function processPlxUnits(datapath)
% Separates merged plx/mwk file into individual files per unit/protocol

    if nargin < 1
        datapath = [pwd, '/'];
    end
    outpath = [datapath, 'offline_data/'];
    if ~exist(outpath, 'dir'); 
        display('Making offline output directory.')
        mkdir(outpath);
    end

    file_d = dir([datapath, '*FULL.mat']);
    for di = 1:length(file_d)
        filename = file_d(di).name;
        filedate_i = union(strfind(filename, '2014'), strfind(filename, '2015'));
        filedate = filename(filedate_i+2:filedate_i+7);
        
        display(filename)
        
        fn = [datapath, filename];
        load(fn, 'RSVP_OBJ', 'RSVP_FBOP', 'M2S');
        if exist('RSVP_OBJ', 'var');
            for ui = 1:length(RSVP_OBJ.allTrials)
                if ~isempty(RSVP_OBJ.allTrials{ui})
                    allTrials = RSVP_OBJ.allTrials{ui};
                    t = findstr(filename, '_');

                    fn2 = ['RSVP', num2str(filedate), '_OBJ_',...
                        filename(t(1)+1:t(2)-1), 'u', num2str(ui-1), '.mat'];
                    
%                     fn2 = ['RSVP', num2str(filedate), '_OBJv2_',...
%                         filename(t(1)+1:t(2)-1), 'u', num2str(ui-1), '.mat'];
                    save([outpath, fn2], 'allTrials');
                    display(['Saved to ', fn2])
                end
            end
        elseif exist('RSVP_FBOP', 'var');
            for ui = 1:length(RSVP_OBJ.allTrials)
                if ~isempty(RSVP_OBJ.allTrials{ui})
                    allTrials = RSVP_OBJ.allTrials{ui}; 
                    t = findstr(filename, '_');

                    fn2 = ['RSVP', num2str(filedate), '_FBOP_',...
                        filename(t(1)+1:t(2)-1), 'u', num2str(ui-1), '.mat'];
                    save([outpath, fn2], 'allTrials');
                    display(['Saved to ', fn2])
                end
            end
            
        elseif exist('M2S', 'var');
            for ui = 1:length(M2S.allTrials)
                if ~isempty(M2S.allTrials{ui})
                    allTrials = M2S.allTrials{ui}; 
                    t = findstr(filename, '_');

                    fn2 = ['M2S', num2str(filedate), '_OBJ_',...
                        filename(t(1)+1:t(2)-1), 'u', num2str(ui-1), '.mat'];
                    save([outpath, fn2], 'allTrials');
                    display(['Saved to ', fn2])
                end
            end
        end
        
        clear RSVP_OBJ RSVP_FBOP M2S
    end

end