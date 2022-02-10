% LLC v3. - Based on analyEVP from MLAnalysisOnline
% January 2022.
% Brock Carlson. 
% 2nd year - Vanderbilt University Psychology Graduate Program

% Master code, from start to finish, of the processing for the Maier/Bastos
% laminar labeling collaboration. This code takes in .ns2 data recorded on
% BlackRock and triggers the LFP to stimulus onsets. Then, we process CSD
% and CSD from the LFP to assign cortical depths from each type of data.
% Finally, we get the average CSD and PSD profile for each penetration,
% aligned to the granular input layer, and compare coherence across
% sessions.

% Desired outputs (per session)
% 1. ch x time x trials -- LFP, PSD, CSD
% 2. ch x time -- LFP, PSD, CSD (with plots of averages)
% 3.  -3000 : 100 : +3000 depth normalized -- PSD CSD.

% final desired output == (penetration x depth x time/freq)

%% FILE SELECTION
clear
close all

%% SET UP DIRECTORIES - I.E. PostSetup('BrockWork')
CODEDIR  = 'C:\Users\Brock\Documents\MATLAB\GitHub\laminarLabelingCollab\MasterScripts';
cd(CODEDIR)
DATADIR = 'E:\all BRFS';
task = 'brfs';
OUTDIR = 'E:\LLC individual penetration outputs';


%% set up directory loop
BRdatafile    = 'E:\all BRFS\151221_E\151221_E_brfs001';

%% Pre-processing the LFP
extension     = 'ns2'; % THIS CODE DOES NOT DOWNSAMPLE OR FILTER DATA
el            = 'eD';
sortdirection = 'ascending'; %  descending (NN) or ascending (Uprobe)
pre           = -200;
post          = 800;

flag_subtractbasline = true;
flag_halfwaverectify = false;


%% Trial-align data -- Create STIM

filelist = {BRdatafile};
V1 = 'LV1'; %% check this later
STIM = diTP(filelist,V1);
[STIM,fails] = diPT(STIM); 

%diV1lim
    global ALIGNDIR
    ALIGNDIR = 'E:\V1Limits\';
    penetrationNumber = 1;
    STIM.penetration = strcat(STIM.header,'_eD');
    STIM.rmch = 0;
    % STIM = diV1Lim(STIM,penetrationNumber);
    TuneList = importTuneList();
    idx = find(strcmp(TuneList.Datestr,STIM.header(1:6)));
    idx = idx(penetrationNumber);
    penetration =  TuneList.Penetration{idx};
    l4c = TuneList.SinkBtm(strcmp(TuneList.Penetration,penetration));
    l4l = sprintf('e%s%02u',penetration(end),l4c);
    aligndir = ALIGNDIR;
    load([aligndir penetration ,'.mat'],'elabel','v1lim','fRF')
    clear v1lim
    v1lim = [1 24];
    el_labels = elabel(v1lim(1):v1lim(2));
    l4_idx    = find(strcmp(el_labels,l4l));
    ninside   = length(el_labels);
    depths    = [0:ninside-1; -1*((1:ninside) - l4_idx); ninside-1:-1:0]';
    STIM.el_labels    =  el_labels;
    STIM.depths       =  depths;
    STIM.v1lim         = [v1lim(1) find(strcmp(elabel,l4l)) v1lim(2)];

% Create LFP triggered SDF
[RESP, win_ms, SDF, sdftm, PSTH, psthtm]= trialAlignLFP_BMC(STIM,pre,post); %note 'pre' must be negative
chans = 1:size(SDF,1);

%% Select trials
SDFch1 = squeeze(SDF(1,:,:));
trls = ~isnan(SDFch1(1001,:));
SDF800 = SDF(:,:,trls);
if sum(trls) < 25
    error('trial count too low. Can you add binocular simultaneous? Or just include trials without a NaN after 500?')
end

% % 
% % EVP = DAT(:,:,trls);
EVP = nanmean(SDF800,3);

%% CSD processing and plotting
switch sortdirection
    case 'ascending' 
        corticaldepth = [chans] ;
    case 'descending'
        corticaldepth = fliplr([chans]);
end

CSD = calcCSD(EVP');
if flag_subtractbasline
    CSD = bsxfun(@minus,CSD,mean(CSD(:,sdftm<0),2));
end
if flag_halfwaverectify
    CSD(CSD > 0) = 0;
end
CSD = padarray(CSD,[1 0],NaN,'replicate');
csdfig = figure;
subplot(1,2,1)
f_ShadedLinePlotbyDepth(CSD,corticaldepth,sdftm,[],1)
title(BRdatafile,'interpreter','none')
set(gcf,'Position',[1 40 700 1200]); 



%%
CSDf = filterCSD(CSD);

subplot(1,2,2)
switch sortdirection
    case 'ascending'
        y = [chans];
        ydir = 'reverse';
    case 'descending'
        y = fliplr([chans]);
        ydir = 'normal';
end
imagesc(sdftm,y,CSDf); colormap(flipud(jet));
climit = max(abs(get(gca,'CLim'))*.8);
set(gca,'CLim',[-climit climit],'Ydir',ydir,'Box','off','TickDir','out')
hold on;
plot([0 0], ylim,'k')
c = colorbar;
% caxis([-250 250])

%% Save figs
saveFigNameFIG = strcat(OUTDIR,filesep,BRdatafile(23:end),'_CSD.fig');
saveFigNamePNG = strcat(OUTDIR,filesep,BRdatafile(23:end),'_CSD.png');
savefig(csdfig,saveFigNameFIG)
saveas(csdfig,saveFigNamePNG)

%% Save data
% avgBeta, avgGamma, power_norm, freq_vector, STIM, powerAvg. SDF800
saveName = strcat(OUTDIR,filesep,BRdatafile(23:end),'_PSD_matVar.m');
save(saveName, 'STIM', 'SDF', 'SDF800', 'CSD', 'sdftm','corticaldepth',...
    '-mat','-v7.3')
