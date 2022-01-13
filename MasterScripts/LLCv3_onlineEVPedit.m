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

%% FILE SELECTION
clear
BRdatafile    = '\\129.59.230.171\CerebrusData\190320_B\190320_B_evp004';

%% Pre-processing
% This section must perform the following tasks:
% 1. load .ns2
% 2. band-stop-filter .ns2 data at 60Hz (skip for now?)
% 3. load .nev
% 4. create "triggerpoints" variable
% 5. create a logical "index" for trials with 800ms of simulation.
% 6. create a trial-triggered variable - "LFP" (ch x time x trls) 

[drname,BRdatafile] = fileparts(BRdatafile); 
filename = fullfile(drname,BRdatafile);
cd(drname);
if exist(strcat(filename,'.nev'),'file') == 2;
    NEV = openNEV(strcat(filename,'.nev'),'nomat','nosave');
else
    error('the following file does not exist\n%s.nev',filename);
end
% get event codes from NEV, then clear
    EventCodes = NEV.Data.SerialDigitalIO.UnparsedData - 128;
    EventTimes = floor(NEV.Data.SerialDigitalIO.TimeStampSec .* 1000); %ms, to match 1kHz
    clear NEV
    
% create labels -- do I need all of this?
NS_Header = openNSx(strcat(filename,'.',extension),'noread');
neural = strcmp({NS_Header.ElectrodesInfo.ConnectorBank},el(2));
N.neural = sum(neural);
NeuralLabels = {NS_Header.ElectrodesInfo(neural).Label};
Fs = NS_Header.MetaTags.SamplingFreq;
% process data electrode by electrode
nct = 0;
for e = 1:length(neural)
    if neural(e) == 1
        % good electrode
        nct = nct+1;
        clear NS DAT
        electrode = sprintf('c:%u',e);
        NS = openNSx(strcat(filename,'.',extension),electrode,'read','uV');
        if iscell(NS.Data)
            DAT = cell2mat(NS.Data);
        else
            DAT = NS.Data;
        end
        NS.Data = [];
        if nct == 1
            %preallocation
            N.samples = length(DAT); %samples in header diffrent from actual data length???
            clear LFP
            LFP_preSortDir = zeros(ceil(N.samples),N.neural);
        end
        LFP_preSortDir(:,nct) = DAT;
        clear DAT
    end
end 

% sort electrode contacts in ascending order:
for ch = 1:length(NeuralLabels)
    chname = strcat(sprintf('%s',el),sprintf('%02d',ch));
    id = find(~cellfun('isempty',strfind(NeuralLabels,chname)));
    if ~isempty(id)
        ids(ch) = id;
    end
end

switch sortdirection
case 'ascending'
    LFP = LFP_preSortDir(:,ids);
case 'descending'
    LFP = LFP_preSortDir(:,fliplr(ids));
    otherwise
    error('need sort direction')
end


%% CSD processing
%   1. calculate CSD
%   2. plot CSD with imagesc
%   3. Identify Layer 4c CSD sink bottom.
%   4. Align all penetrations along depth via CSD sink bottom.


%% PSD processing
%   1. calculate PSD
%   2. plot PSD with imagesc and gamma x beta cross
%   2. visually identify gamma x beta cross.
%   3. Align all penetrations along depth via y x b cross.

