clear
% close all
BRdatafile    = 'D:\all BRFS\151231_E\151231_E_brfs001'
extension     = 'ns2'; % THIS CODE DOES NOT DOWNSAMPLE OR FILTER DATA
el            = 'eD';
sortdirection = 'ascending'; %  descending (NN) or ascending (Uprobe)
pre           = 200;
post          = 800;
chans         = [1:24]; 
trls          = [1:200];


flag_subtractbasline = true;
flag_halfwaverectify = false;

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

[DAT, TM] = trigData(LFP, triggerpoints , pre, post);
if isempty(trls)
    EVP = mean(DAT,3);
else
EVP = mean(DAT(:,:,trls),3);
end

        
%%
%figure;
switch sortdirection
    case 'ascending' 
        corticaldepth = [chans] ;
    case 'descending'
        corticaldepth = fliplr([chans]);
end
%f_ShadedLinePlotbyDepth(EVP,corticaldepth,TM,[],1)
%title(BRdatafile,'interpreter','none')

%%
CSD = calcCSD(EVP);
if flag_subtractbasline
    CSD = bsxfun(@minus,CSD,mean(CSD(:,TM<0),2));
end
if flag_halfwaverectify
    CSD(CSD > 0) = 0;
end
CSD = padarray(CSD,[1 0],NaN,'replicate');
figure
subplot(1,2,1)
f_ShadedLinePlotbyDepth(CSD,corticaldepth,TM,[],1)
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
imagesc(TM,y,CSDf); colormap(flipud(jet));
climit = max(abs(get(gca,'CLim'))*.8);
set(gca,'CLim',[-climit climit],'Ydir',ydir,'Box','off','TickDir','out')
hold on;
plot([0 0], ylim,'k')
c = colorbar;
% caxis([-250 250])

% %%
% clear
% addpath(genpath('C:\Users\maierlab\Documents\MLAnalysisOnline')); 
% BRdatafile    = '\\129.59.230.171\CerebrusData\190313_B\190313_B_evp002';
% extension     = 'ns2'; % THIS CODE DOES NOT DOWNSAMPLE OR FILTER DATA
% sortdirection = 'ascending'; 
% el            = 'eD';
% signals          = {'mua'}; 
% theseelectrodes = {'eD02';'eD03';'eD04';'eD05';'eD06';'eD07';'eD08';'eD09';'eD10';'eD11'; 'eD12';'eD13';'eD14';'eD15';'eD16';'eD17';'eD18';'eD19';'eD20';'eD21';'eD22';'eD23';'eD24','eD25','eD26','eD27','eD28','eD29','eD30','eD31'}; 
% [~, EventCodes, EventTimes]= getLFP(BRdatafile,extension,el,sortdirection);
% [LFP,MUA] = getNeuralData(signals,theseelectrodes,BRdatafile); 
% triggerpoints = EventTimes(EventCodes == 23 | EventCodes == 25 | EventCodes == 27 | EventCodes == 29| EventCodes == 31);
% 
% muatr = trigData(MUA,floor(triggerpoints./30),100,200); 
% 
% figure, plot(-100:200,nanmean(muatr,3)); 