function LFP_bsf = bandStopFiltLFP(LFP)


Fs = 1000;
nyq = Fs/2;

bandStop = [59 61];  %band stop values in Hz
bandStopRad = bandStop*2*pi./Fs; %convert to rad/sample
          %rad/sample = cycles/sec * 1sec/1000samples * radians(2pi)/cycle
[bwb,bwa] = butter(4,bandStopRad,'stop');
    % want to visualize this filter? Use --> freqz(bwb,bwa)
stopLFP = filtfilt(bwb,bwa,LFP); %bandStopFilter60Hz


fullRectLFP = abs(stopLFP);


LFP_bsf = fullRectLFP;


