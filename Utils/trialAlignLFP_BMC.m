function [RESP, win_ms, SDF, sdftm, PSTH, psthtm, SUA, spktm]= trialAlignLFP_BMC(STIM,pre,post)

    datatype = 'lfp';
    win_ms    = [50 100; 150 250; 50 250; -50 0]; % ms;

    flag_1kHz = true;




%% 1. Setup for this penetration   


tfdouble = isa(STIM.tp_pt,'double');
if ~tfdouble
    warning('Expecting STIM.tp_pt to be a double but it is not. Unexpected behavior possible.')
end

global AUTODIR
if ~isempty(AUTODIR)
    autodir = AUTODIR;
else
    autodir = '/Volumes/Drobo2/DATA/NEUROPHYS/AutoSort-ed/';
end
global SORTDIR
if ~isempty(SORTDIR)
    sortdir = SORTDIR;
else
    sortdir = '/Volumes/Drobo2/DATA/NEUROPHYS/KiloSort-ed/';
end

bin_ms    = 10;
nwin      = size(win_ms,1);

el_labels = STIM.el_labels;

nel = length(el_labels);



%% 2.  Preallocate nans for all datatypes (sdftm and triggering is setup)

Fs = 1000;
TP = round(STIM.tp_pt ./ 30); %DEV: make TP/r relationship clearer, it's handeled diffrently for spk data
r = 1; 

RESP      = nan(nel,nwin,length(STIM.trl));

sdftm     = [-0.3*Fs: 0.15*Fs + max(diff(TP,[],2))]; 
SDF       = nan(nel,length(sdftm),length(STIM.trl));

psthtm     = [];
PSTH      = [];
    

%% 3. Load in data, overwriting NaNs to all output formats (SDF, PSTH, etc)
% Iterate trials, loading files as you go
for i = 1:length(STIM.trl)
    
    if i == 1 || STIM.filen(i) ~= filen 
        
        clear filen filename BRdatafile
        filen = STIM.filen(i);
        filename  = STIM.filelist{filen};
        [~,BRdatafile,~] = fileparts(filename);
        
        % setup SPK cell array
        SPK   = cell(nel,1);
        empty = false(size(SPK)); 
        
           
        clear ns2file ns_header
        ns2file    = [filename '.ns2'];
        ns_header  = openNSx(ns2file,'noread');
        
        clear elb idx e
        elb = cellfun(@(x) (x(1:4)),{ns_header.ElectrodesInfo.Label},'UniformOutput',0);
        idx = zeros(1,length(el_labels));
        for e = 1:length(el_labels)
            idx(e) = find(strcmp(elb,el_labels{e}));
        end 
    end
    
   
    
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
    

     for w = 1:nwin
         if rwin(w,1) == rwin(w,2)
             continue
         end
         timeperiod = sprintf('t:%u:%u', rwin(w,1),rwin(w,2));
         NS = openNSx(ns2file,timeperiod,...
             'read','sample');
         dat = double((NS.Data(idx,:)))' ./ 4;  clear NS; 
         RESP(:,w,i) = nanmean(dat,1);
     end
    
     if tp(1) + sdftm(end) > ns_header.MetaTags.DataPoints 
         continue
     end
     timeperiod = sprintf('t:%u:%u',tp(1) + [sdftm(1) sdftm(end)]);
     NS = openNSx(ns2file,timeperiod,...
         'read','sample');
     dat = double((NS.Data(idx,:))) ./ 4; clear NS; 
     dat(:,sdftm > diff(tp)) = [];
     SDF(:,1:length(dat),i) = dat;

      
            
end % done iterating trials

%% 4. Cut/trim so that everything lines up.
% remove last bin (inf) and center time vector
if ~isempty(PSTH)
    PSTH(:,end,:) = [];
    psthtm(end) = [];
    psthtm = psthtm + bin_sp/2;
    % convert time vector to seconds
    psthtm = psthtm./Fs;
end




% trim SDF of convolution extreams %MAY NEED DEV
trim = sdftm < pre | sdftm > post;
sdftm(trim) = [];
SDF(:,trim,:) = [];

% convert time vector to seconds
sdftm = sdftm./(Fs/r);


% trim PSTH to match
if ~isempty(PSTH)
    trim = psthtm < -0.15 | psthtm > psthtm(end) -0.15;
    psthtm(trim) = [];
    PSTH(:,trim,:) = [];
end
