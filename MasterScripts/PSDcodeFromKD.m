
clear
% close all
BRdatafile    = 'D:\all BRFS\151231_E\151231_E_brfs001'
extension     = 'ns2'; % THIS CODE DOES NOT DOWNSAMPLE OR FILTER DATA
el            = 'eD';
sortdirection = 'ascending'; %  descending (NN) or ascending (Uprobe)
chans         = 1:24; 

[LFP, EventCodes, EventTimes]= getLFP(BRdatafile,extension,el,sortdirection);
triggerpoints = EventTimes(EventCodes == 23 | EventCodes == 25 | EventCodes == 27 | EventCodes == 29| EventCodes == 31);

lfp = LFP(:,chans);


Fs       = 1000; % Hz
chanN    = size(lfp,2); 

% loop through channels 
for ch = 1:chanN
clear x n Spec
x        = lfp(:,ch);   
n        = size(lfp,1); % Number of data points
% prep for psd 
    nfft     = 512; 
    window   = hanning(nfft); 
    nwind    = length(window); 
if n < nwind    % zero-pad x if it has length less than the window length
    x(nwind)=0;  n=nwind;
end
noverlap = 1;
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

% normalize power @ each frequency relative to power across contacts 
 power_norm = nan(size(power)); 
 for ch = 1:size(power,2)
     for f = 1:size(power,1)
         power_norm(f,ch) = (power(f,ch) - mean(power(f,:)))./(mean(power(f,:))) * 100; % percent deviation from mean
     end
 end

figure(1), set(gcf,'color','w','position',[1 1 400 800]); 
imagesc(freq_vector,1:24,power_norm'); 
colormap('hot'); xlim([0 100]); 
xlabel('freq (Hz)'); ylabel('contact number'); 
set(gca,'tickdir','out','ytick',1:24); 

%% Get the Gamma x Beta cross
% Beta is 12 - 20Hz (for our purposes)
% Gamma is 30-59,61:100
figure(2)
beta_index = (freq_vector > 12) & (freq_vector < 25);
gamma_index = (freq_vector > 30);

clear i avgBeta avgGamma
for i = chans
    avgBeta(i,1) = mean(power_norm(beta_index,i));
    avgGamma(i,1) = mean(power_norm(gamma_index,i));
end


plot(avgBeta)
hold on
plot(fliplr(avgGamma))
view([90 -90])
set(gca,'xdir','reverse')
legend('Beta','Gamma','Location','best')
titleText = {'Normalized Gamma x Beta power across contacts',BRdatafile(22:end)};
title(titleText,'Interpreter','none')