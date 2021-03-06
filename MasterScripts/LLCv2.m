% LLC v2.
% December 2021.
% Brock Carlson. 
% 2nd year - Vanderbilt University Psychology Graduate Program

% Master code, from start to finish, of the processing for the Maier/Bastos
% laminar labeling collaboration. This code takes in .ns2 data recorded on
% BlackRock and triggers the LFP to stimulus onsets. Then, we process CSD
% and CSD from the LFP to assign cortical depths from each type of data.
% Finally, we get the average CSD and PSD profile for each penetration,
% aligned to the granular input layer, and compare coherence across
% sessions.

%% Pre-processing
% This section must perform the following tasks:

%   1. Process the .NEV file to find relevant stimulus onsets.
%       Relevant stimuli are any stimuli on the screen for at least 500ms. 
%       The macaque subject must have fixated through the whole trial.
%   2. Align the .ns2 (1kHz) LFP data to stimulus onsets
%       We then subsequently align to the photo diode onset

clear
setup_llc('home')

% runTuneList.m breakdown
TuneList = importTuneList(3);
paradigm = {'cosinteroc','mcosinteroc','brfs'};

% cut down TuneList to only analyze files from a given RIGDIR if needed for testing.    
if flag_UseRigDir
    global RIGDIR
    list = dir([RIGDIR]);
    dirlist = cellfun(@(x,y) sprintf('%s_%s',x,y),TuneList.Datestr,TuneList.Monkey,'UniformOutput',false);

    I = ismember(dirlist,{list.name}); % logical output
     fields = fieldnames(TuneList);
        for f = 1:length(fields)
            TuneList.(fields{f})(~I) = [];
        end
end




% ADD HERE!! bp filter the ns2 file, save as ns2 to new file?
for s = 1:length(TuneList.Penetration) 

    %filter dat here
        % get TP and rwin
    clear tp rwin
    tp = TP(i,:) ; % TP is always Fs = 30000 because it is recorded by NEV
    if any(isnan(tp))
        continue
    end
    
    rwin = tp(1) + (win_ms/1000*Fs);
    % deal with instnaces where the window exceeds stimulus offset
    stimoff = rwin > tp(2); 
    rwin(stimoff)  = tp(2); 
    
    
    
        NS = openNSx(filename,'read','sample');
    dat = double(NS.Data) ./ 4; clear NS; 
    dat(:,sdftm > diff(tp)) = [];
    
    
    save dat. -- can I recover the timeponits?

end


    % example
    filename = 'E:\all BRFS\151221_E\151221_E_brfs001';
    ns2file    = [filename '.ns2'];
    tic
    NS = openNSx(ns2file,'read','sample'); %in 1kHz samples. I can use same TP.
    toc
    dat = double(NS.Data) ./ 4; %weird data conversion needed. Converts to uV?
    clear NS; 




% Get timepoints
for s = 1:length(TuneList.Penetration) 
    add file sorting stuff for penetration here
    
    % get event codes from NEV
        STIM = diTP(filelist,V1);
    
    % get photodiode signal from .ns6 file ??Can I use ns2 file here?
        add some diPT suff here (to cut the diPT file)
        [newTP,message] = photoReTriggerSTIM(...
                                STIM.tp_sp(I,:),...
                                filename,...
                                STIM.ypos(I,:));
                            
    
    % Set sink bottom (and cut cortex?) (can I include all channels here?)
    	STIM = diV1Lim(STIM,pn);


    global STIMDIR    
    save([STIMDIR STIM.penetration '.mat'],...
        'STIM', '-v7.3')
    clear STIM

        
    
end





%%

% diRunDir breakdown
global STIMDIR
cd(didir)
list = dir([didir '2*.mat']); % need a better way to Identify .mat files we need. 
I = cellfun(@length,{list.name}) == 15;
list = list(I);  %creates a list of STIM.mat files 
for i = 1:length(list)

    load([didir list(i).name],'STIM')
    datastr = upper(datatype);
    matname = [didir list(i).name(1:end-4) '_' datastr '.mat']; 
    [RESP, win_ms, SDF, sdftm, PSTH, psthtm] = diNeuralDat(STIM,datatype,true);
            within this I shoudl bp filter all data. No timeperoid.
% %              for w = 1:nwin
% %                  if rwin(w,1) == rwin(w,2)
% %                      continue
% %                  end
% %                  timeperiod = sprintf('t:%u:%u', rwin(w,1),rwin(w,2));
% %                  NS = openNSx(ns2file,timeperiod,...
% %                      'read','sample');
% %                  dat = double((NS.Data(idx,:)))' ./ 4;  clear NS; 
% %                  RESP(:,w,i) = nanmean(dat,1);
% %              end
    
    save(matname,...
        'STIM','RESP','win_ms','SDF','sdftm','PSTH','psthtm',...
        '-v7.3')
    
end

%%%% GOAL OUTPUT: anaType = '_LFP.mat'; in STIMDIR




%% IDX processing (laminarLabeling_LFP.m)
% Load the penetration_STIM_LFP.mat files
% Loop through filelist
    % Loop through each electrode
        % Pull out SDF for trials that have data out to 500ms.
        % crop/pad SDF
        % Save IDX


%%  Reformat LFP to fix grand matrix on an individual session level (reformatLamLabLFP.m)
%%%%% YOU CAN INCLUDE THIS IN THE PREVIOUS SECTION AND MAKE IT FASTER?
% format the LFP (ch x time x trls) for each session
% export


%% CSD processing
%   1. Identify Layer 4c CSD sink bottom.
%   2. Align all penetrations along depth via CSD sink bottom.


%% PSD processing
%   1. Identify gamma x beta cross
%   2. Align all penetrations along depth via PSD's gamma x beta cross

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Pre-processing - utilizing the ephys analysis framework.
% This section must perform the following tasks:

%   1. Process the .NEV file to find relevant stimulus onsets.
%       Relevant stimuli are any stimuli on the screen for at least 500ms. 
%       The macaque subject must have fixated through the whole trial.
%   2. Align the .ns2 (1kHz) LFP data to stimulus onsets
%       We then subsequently align to the photo diode onset



% Recreate analysis done in runTuneList. Do not crop v1lim.
% 1. importTuneList - establish 


% 2. loop through penetrations
    % 2a. pull info from TuneList
    % 2b. get TPs
    % 2c. diCheck (can we skip this?)
    % 2d. photodiode trigger
    % 2e. skip diV1Lim but do establish the STIM.el_labels and STIM.rank
    
    
%%%%% Can I skip all of this and use the 1.75 GB of LFP stim on my PC? It
%%%%% is not cut by cortical depth but we can go back and grab that later
%%%%% if we need. 


%% Get properly formatted LFP (ch x time x trial) for each session
createSTIM_LFP = false;

clear
close all
PostSetup('BrockHome_LLC')

if createSTIM_LFP
    global IDXDIR
    cd(IDXDIR)
    if ~exist(strcat(IDXDIR,'\laminarLabeling_LFP.mat'),'file')
        laminarLabeling_LFP
    end
        load(strcat(IDXDIR,'\laminarLabeling_LFP.mat')) %317s (aka 5 min)
    % Save the individual LFP session-wise matrices
        for i = 1:size(IDX.allV1,2)
            sessions{i} = IDX.allV1(i).penetration; 
        end
        individualSes = unique(sessions)';
        clearvars -except IDX
    % 30 sessions. 25 from E, 5 from I
        lfpDir = 'E:\5 diIDX dir\laminarLabelingLFPs\';
        cd(lfpDir)
        refLLL = false;
        if refLLL
            reformatLamLabLFP(IDX) %138sec (aka 3.3min)
        else
            disp('you have LFP properly formatted - no need to run')
        end
        clear IDX
end

% Examine LFP outputs and make sure they look reasonable. All should be bl
% corrected.
global LFPDIR
cd(LFPDIR)
anaType = '_LFP_4LLC.mat';
list    = dir([LFPDIR '*' anaType]);
close all
for i = 1:length(list)
    clear penetration
    penetration = list(i).name; 
    clear LFP avgLFP
    load(penetration)
    % Average LFP at each contact
        avgLFP = mean(LFP,3)'; % we flip the first two dimension outputs to plot
    % plot
        TM = -.15:.001:.499;
        figure
        plot(TM,avgLFP); hold on; 
        vline(0)
        hline(0)
        xlim([-.15 .5])
        xlabel('time (in ms)')
        ylabel('uV')
        title(penetration,'interpreter','none')
        
        
    
end
%% Calculate CSD from the LFP on a session level
% 1. replicate top and bottom lfp channels
% 1.5 Calculate CSD.
% 2. Z-score normalize each CSD channel across trials
% 3. Average across trials. Output is (ch x time).
% 4. Format depth in terms of grand alignment depth (100um).
% 5. Save each session's CSD.

%% Create CSD master matrix
% Pull out all sessions of CSD data and concatenate along 1st dimension.
% output is (session x ch x time)

%24 contqacts at spacing of 100um
coordinate_axis = -3200:100:3200;
dist_to_top= sink-1;
dist_to_bottom = 24-sink+1;
top_coord = dist_to_top*-100;
bottom_coord = dist_to_bottom*100;
new_coords = top_coord:100:bottom_coord; %cop into ElectrodeInfo
indx1 = find(coordinate_axis== top_coord)
indx2 = find(coordinate_axis== bottom_coord)
master_matrix_POW(n,indx1:indx2,:) = this_session_POW;
master_matrix_CSD(n,indx1:indx2,:) = this_session_CSD;



