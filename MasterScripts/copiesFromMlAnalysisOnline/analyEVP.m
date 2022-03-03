clearvars -except ASubBase
BRdatafile    = 'T:\rig021_LaminarLabelingCollaboration\EndOfDayFileOutputs\220228_B\220228_B_evp002';

extension     = 'ns2'; % THIS CODE DOES NOT DOWNSAMPLE OR FILTER DATA
el            = 'eA';
sortdirection = 'ascending'; %  descending (NN) or ascending (Uprobe) % new note -- BMC 211007_B descenting and ascending is a moot point because of the new channel map
% pre           = 50;
% post          = 250;
pre           = 50;
post          = 250;
%chans         = [1:17 19:20 22:32];
chans    =    [1:32]; 
%chans           = [1:2:24];
trls          = [1:1000];

 
flag_subtractbasline = true;
flag_halfwaverectify = false;
flag_normalize = false;
flag_interpolate = true;
interp_chans = [5,12,19];

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

if flag_normalize
    EVP = normalize(EVP,2);
end

if flag_interpolate 
    for i = 1:length(interp_chans)
        badChan = interp_chans(i);
        EVP(:,badChan) = (EVP(:,badChan+1) + EVP(:,badChan-1)) / 2;
    end
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
set(gcf,'Position',[1 40  700   785]); 
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
caxis([-500 500])

