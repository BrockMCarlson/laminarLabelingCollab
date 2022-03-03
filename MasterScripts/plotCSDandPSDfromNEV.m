%% Plot CSD and PSD with NEV 
% no .bhv file required
% Combined from analyEVP.m and analyPSDdepth_redo.m in MLAnalysisOnline


clear
close all

%% Set up file list
% dataFileLocation = 'T:\rig021_LaminarLabelingCollaboration\EndOfDayFileOutputs\';
dataFileLocation = '\\CEREBUSHOSTPC\CerebrusData\';
folderName = ['220131_B';"220202_B";"220204_B";"220207_B";"220209_B";...
    "220211_B";"220214_B";"220216_B";"220218_B";"220221_B";"220223_B";...
    "220225_B";"220228_B"];
evpNumber = ['4';'3';'5';'7';'3';'5';'3';'8';'3';'4';'4';'4';'2'];
for i = 1:13
    fullFileName(i,1) = ...
        strcat(dataFileLocation, folderName(i), filesep, folderName(i),...
        '_evp00', evpNumber(i));
end
useChans = {1:24; 1:24; 1:24; 1:31; 1:32; 1:32; 1:32; 1:32; 1:32; 1:32; 1:32; 1:32; 1:32};
interpTheseChans = {[15,22]; [15,22]; [15,22]; [14]; [18]; [5,10]; [16]; [13 16]; [16]; [16]; []; []; []};
useSession = [false; false; false; false; true; true; true; true; true; true; true; true; true];

FileInformation = table(folderName,useSession,evpNumber,useChans,interpTheseChans,fullFileName);
clearvars -except FileInformation

%% Choose your session number
% SessionNum = 8;
for SessionNum = 1:13
%%
if ~FileInformation.useSession(SessionNum)
    continue
end

BRdatafile = FileInformation.fullFileName{SessionNum};
extension     = 'ns2'; % THIS CODE DOES NOT DOWNSAMPLE OR FILTER DATA
el            = 'eA';
sortdirection = 'ascending'; %  descending (NN) or ascending (Uprobe) % new note -- BMC 211007_B descenting and ascending is a moot point because of the new channel map
pre           = 50;
post          = 250;
chans = FileInformation.useChans{SessionNum};

 
flag_subtractbasline = true;
flag_halfwaverectify = false;
flag_interpolate = true;
interp_chans = FileInformation.interpTheseChans{SessionNum};

clear LFP EventCodes EventTimes DAT TM CSD CSDf corticaldepth y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[LFP, EventCodes, EventTimes]= getLFP(BRdatafile,extension,el,sortdirection);
triggerpoints = EventTimes(EventCodes == 23 | EventCodes == 25 | EventCodes == 27 | EventCodes == 29| EventCodes == 31);
if isempty(chans)
    chans = [1:size(LFP,2)];
end
LFP = LFP(:,chans);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[DAT, TM] = trigData(LFP, triggerpoints , pre, post); %DAT comes out in form of (time x channels x trials
DATpermute = permute(DAT,[2 1 3]); %dim chould be chan x ms

EVP = mean(DATpermute,3);


if flag_interpolate 
    for i = 1:length(interp_chans)
        badChan = interp_chans(i);
        EVP(badChan,:) = (EVP(badChan+1,:) + EVP(badChan-1,:)) / 2;
    end
end


        
%%
%figure;
switch sortdirection
    case 'ascending' 
        corticaldepth = chans ;
    case 'descending'
        corticaldepth = fliplr(chans);
end
%f_ShadedLinePlotbyDepth(EVP,corticaldepth,TM,[],1)
%title(BRdatafile,'interpreter','none')

%%
CSD = calcCSD_classic(EVP); % contains padArray. Input is time x chan. Output is Chan x time
if flag_subtractbasline
    CSD = bsxfun(@minus,CSD,mean(CSD(:,TM<0),2));
end
if flag_halfwaverectify
    CSD(CSD > 0) = 0;
end
figure
subplot(1,5,1)
f_ShadedLinePlotbyDepth(CSD,corticaldepth,TM,[],1)
sgtitle(BRdatafile,'interpreter','none')
%%
CSDf = filterCSD(CSD);

subplot(1,5,2)
switch sortdirection
    case 'ascending'
        y = [chans];
        ydir = 'reverse';
    case 'descending'
        y = fliplr([chans]);
        ydir = 'normal';
end
imagesc(TM,y,CSDf); colormap(flipud(jet));
climit = max(abs(get(gca,'CLim'))*.8);
set(gca,'CLim',[-climit climit],'Ydir',ydir,'Box','off','TickDir','out')
hold on;
plot([0 0], ylim,'k')
c = colorbar;
caxis([-500 500])


%% PSD code

%%
%%


chanLim = chans;
jnmfile = [BRdatafile '.ns2'];
badchan = [];

if ~exist( 'val', 'var' )    
    val = 15; 
end

t = openNSx( jnmfile, 'noread' );

tn = ~strcmp( 'E', { t.ElectrodesInfo.ConnectorBank } );
nsx.E = length( tn );
nsx.neural = sum( tn );
nsx.anlg = sum( ~tn );

if exist('electrode')
error('what goes here?')   
end

neuralChan = find( tn );
bncChan = find( ~tn );

nsx.neuralL = { t.ElectrodesInfo( tn ).Label };
nsx.neuralI = t.ElectrodesInfo( tn );

nsx.bncL = { t.ElectrodesInfo( ~tn ).Label };
nsx.bncI = t.ElectrodesInfo( ~tn );

nsx.rnge = ...
    ( ( length( t.ElectrodesInfo( 1 ).MinAnalogValue : ...
    t.ElectrodesInfo( 1 ).MaxAnalogValue ) ) ./ ...
    ( length( t.ElectrodesInfo( 1 ).MinDigiValue : ...
    t.ElectrodesInfo( 1 ).MaxDigiValue ) ) );

nsx.fs = t.MetaTags.SamplingFreq;
nsx.nyq = nsx.fs / 2;
nsx.deci = nsx.fs / 1000;

electD = openNSx( jnmfile, 'c:1', 'read' );
tData = double( electD.Data );

samples = length( tData );

bnc = zeros( nsx.anlg, ceil( samples / nsx.deci )  );
lfp = zeros( nsx.neural, ceil( samples / nsx.deci ) );
    
electD = openNSx( jnmfile );

tData = ( double( electD.Data ) ).' ;
t2Data = tData( :, bncChan );
tData = tData( :, neuralChan );

electD.Data = [];

tData = tData ./ 4;

lfp = tData.';
bnc = t2Data.';

clear tData t2Data electD 

if size( lfp, 1 ) < 33
    
    pNum = 1;
    
else
    
    pNum = 2;
    
end


if pNum == 2
    pChan(1) = 24;
    idx(1,:) = [1 24];
    pChan(2) = 32;
    idx(2,:) = [25 56];
else
    
pChan = size( lfp, 1 ) / pNum;
idx(1,:) = [1 pChan];
end


    

for k = 1 : pNum
    
lfp2 = jnm_reorder( lfp(idx(k,1):idx(k,2),:), nsx.neuralL(idx(k,1):idx(k,2)), 'BR', pNum, pChan(k) );

if flag_interpolate 
    for i = 1:length(interp_chans)
        badChan = interp_chans(i);
    lfp2(badChan,:) = (lfp2(badChan+1,:) + lfp2(badChan-1,:)) / 2;
    end
end


valn = 2^val;

    jnm = zeros( numel(chanLim), 257 );
    ictr = 1;
    for i = chanLim
        [ jnm( ictr, : ),~,jnmf(ictr,:) ] = jnm_psd( lfp2( i, end - valn : end ), ...
            512, 1000, 512, 0);
        ictr = ictr+1;
    end
    
    jnm2 = zeros( numel(chanLim), 52 );
    
    for i = 1 : 52
        
        jctr = 1;
        for j = chanLim
            
            jnm2( jctr, i ) = ( jnm( jctr, i ) - mean( jnm( :, i ) ) ) ...
                / mean( jnm( :, i ) ) * 100;
            jctr =jctr+1;
        end
    end
    
    %jnm2(:,31:33) = 0; %Manually remove 60Hz signal
    
    
end


[~,fname] = fileparts(BRdatafile);
title(fname,'interpreter','none')

%frequency axis is jnmf(77) = 148Hz




jnmrm = jnm;
jnmrm(badchan,:) = NaN;

normpowAB = jnmrm(:,:) ./ nanmax(jnmrm(:,:), [], 1);

subplot(1,5,[3 4]);
imagesc(normpowAB(:,1:100))

psdAx = gca;
colormap(psdAx,'Jet')

xlabel ={};
indx = 1;
for f = 1:5:100
    tmp = jnmf(1,f);
    tmp = round(tmp);
    xlabel{indx} = num2str(tmp);
    indx=indx+1;
end

set(gca, 'xtick', 1:5:100);
set(gca, 'xticklabel', xlabel)
d = colorbar;
set(gcf,'Position',[5.8000 98.6000 1.5072e+03 671.2000]); 

xlim([1 100])

%% Gamma beta cross
freq_vector = jnmf(1,:);
beta_index = (freq_vector > 10) & (freq_vector < 35);
gamma_index = (freq_vector > 38) ;

clear i avgBeta avgGamma
for i = 1:size(normpowAB,1)
    avgBeta(i,1) = mean(normpowAB(i,beta_index));
    avgGamma(i,1) = mean(normpowAB(i,gamma_index));
end

subplot(1,5,5)
plot(avgBeta)
hold on
plot(fliplr(avgGamma))
view([90 -90])
set(gca,'xdir','reverse')
legend('Beta','Gamma','Location','best')
titleText = {'Normalized Gamma x Beta power across contacts',BRdatafile(22:end)};


%% Save variables
allCSD(SessionNum,:,:) = CSD;
allPSD(SessionNum,:,:) = normpowAB;

end