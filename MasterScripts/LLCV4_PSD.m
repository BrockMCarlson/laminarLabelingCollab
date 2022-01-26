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
BRdatafile    = 'E:\all BRFS\151221_E\151221_E_brfs001';

%% Pre-processing the LFP
extension     = 'ns2'; % THIS CODE DOES NOT DOWNSAMPLE OR FILTER DATA
el            = 'eD';
sortdirection = 'ascending'; %  descending (NN) or ascending (Uprobe)
pre           = 200;
post          = 800;

flag_subtractbasline = true;
flag_halfwaverectify = false;

%% Trigger data
[lfp, EventCodes, EventTimes]= getLFP(BRdatafile,extension,el,sortdirection);
triggerpoints = EventTimes(EventCodes == 23 | EventCodes == 25 | EventCodes == 27 | EventCodes == 29| EventCodes == 31);


[DAT, TM] = trigData(lfp, triggerpoints , pre, post);

%% Select trials
% % trls % you need to pull out the right trials here
% Import .txt file here.

% % EVP = DAT(:,:,trls);


%% FFT
Fs       = 1000; % Hz
chanN    = size(lfp,2); 
% loop through channels 
for ch = 1:chanN
    clear x n Spec
    lfp_holder        = lfp(:,ch);

    %%%%%%%%%%%%%%%%%
    % % % % % % % % % %  x =bandStopFiltLFP(lfp_holder); fix bandstop filter
    %%%%%%%%%%%%%%%%%

    n        = size(lfp,1); % Number of data points
    % prep for psd 
        nfft     = 512; 
        window   = hanning(nfft); 
        nwind    = length(window); 
    if n < nwind    % zero-pad x if it has length less than the window length
        x(nwind)=0;  n=nwind;
    end
    noverlap = 0;
    k        = fix((n-noverlap)/(nwind-noverlap));	% Number of windows
    index    = 1:nwind;
    % compute PSD
        x        = lfp(:,ch);   
        Spec     = zeros(nfft,1); 
    % loop through windows 
        for i=1:k
            xw    = window.*(x(index));
            index = index + (nwind - noverlap);
            Xx    = abs(fft(xw,nfft)).^2;
            Spec  = Spec + Xx;  
        end
    % Select first half
        if ~any(any(imag(x)~=0))   % check if x is complex 
            if rem(nfft,2)    % nfft odd
                select = (1:(nfft+1)/2)';
            else
                select = (1:nfft/2+1)';
            end
            Spec = Spec(select);
        else
            select = (1:nfft)';
        end
    freq_vector = (select - 1)*Fs/nfft;
    if ch == 1
        power = nan(size(Spec,1),chanN); 
    end
    power(:,ch) = Spec; 

end

%% cheat at your band stop filter
idx60hz = find((freq_vector > 57.4 & freq_vector < 62.6 ));
power(idx60hz,:) = 0;

% normalize power @ each frequency relative to power across contacts 
 power_norm = nan(size(power)); 
 for ch = 1:size(power,2)
     for f = 1:size(power,1)
         power_norm(f,ch) = (power(f,ch) - mean(power(f,:)))./(mean(power(f,:))) * 100; % percent deviation from mean
     end
 end

 chans = 1:size(power_norm,2);
figure, set(gcf,'color','w','position',[1 1 400 800]); 
imagesc(freq_vector,chans,power_norm'); 
colormap('hot'); xlim([0 100]); 
xlabel('freq (Hz)'); ylabel('contact number'); 
set(gca,'tickdir','out','ytick',chans); 

%% Get the Gamma x Beta cross
% Beta is 12 - 20Hz (for our purposes)
% Gamma is 30-59,61:100
figure
beta_index = (freq_vector > 12) & (freq_vector < 25);
gamma_index = (freq_vector > 30) ;
gamma_index(idx60hz) = false;

clear i avgBeta avgGamma
for i = chans
    avgBeta(i,1) = mean(power_norm(beta_index,i));
    avgGamma(i,1) = nanmean(power_norm(gamma_index,i));
end


plot(avgBeta)
hold on
plot(fliplr(avgGamma))
view([90 -90])
set(gca,'xdir','reverse')
legend('Beta','Gamma','Location','best')
titleText = {'Normalized Gamma x Beta power across contacts',BRdatafile(22:end)};
