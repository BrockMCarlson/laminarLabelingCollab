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
OUTDIR = 'E:\LLC individual penetration outputs';
ALIGNDIR = 'E:\V1Limits\';

%% Pre-processing the LFP
extension     = 'ns2'; % THIS CODE DOES NOT DOWNSAMPLE OR FILTER DATA
el            = 'eD';
sortdirection = 'ascending'; %  descending (NN) or ascending (Uprobe)
pre           = -200;
post          = 800;

flag_subtractbasline = true;
flag_halfwaverectify = false;

%% set up directory loop
list = dir(DATADIR);
filelist = list(3:end);

for h = 1:size(filelist,1)
    clearvars -except CODEDIR DATADIR OUTDIR ALIGNDIR...
        extension el sortdirection pre post...
        flag_subtractbasline flag_halfwaverectify...
        h filelist
    
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SECTION 0 -- loop through brfs recordings for day %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % BRdatafile
        containingFolder = filelist(h).folder;
        recordingDateName = filelist(h).name;
        recordingDateFullFolder = strcat(containingFolder,filesep,recordingDateName);
        sublist = dir(recordingDateFullFolder);
        brfsIdx = find(...
                       (contains({sublist.name}, 'brfs')) & ...
                       (contains({sublist.name}, 'ns2')) ...
                                                            ==1) ;
        for fileLoop = 1:size(brfsIdx,2)
            ns2BRFSfiles{fileLoop} =...
                strcat(recordingDateFullFolder,filesep,sublist(brfsIdx(fileLoop)).name(1:end-4));
        end
    
    
    for i = 1:size(ns2BRFSfiles,2)
        try
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SECTION 1 -- Create STIM -- Load LFP and Trial-Align %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        close all

        % BRdatafile
            BRdatafile = ns2BRFSfiles{i};
        
        % diTP - create STIM with condition information
            V1 = 'LV1'; %% check this later
            STIM = diTP({BRdatafile},V1);
            [STIM,fails] = diPT(STIM); 
        
        %diV1lim
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
            load([ALIGNDIR penetration ,'.mat'],'elabel','v1lim','fRF')
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
        
        % Select trials
            SDFch1 = squeeze(SDF(1,:,:));
            trls = ~isnan(SDFch1(1001,:));
            SDF800 = SDF(:,:,trls);
            if sum(trls) < 25
                error('trial count too low. Can you add binocular simultaneous? Or just include trials without a NaN after 500?')
            end
            
        
            EVP = mean(SDF800,3);
        
        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SECTION 2 -- CSD processing and plotting %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        % 
        switch sortdirection
            case 'ascending' 
                corticaldepth = chans ;
            case 'descending'
                corticaldepth = fliplr(chans);
        end
        
        %CalcCSD -- 
        CSD = calcCSD(EVP');
        if flag_subtractbasline
            CSD = bsxfun(@minus,CSD,mean(CSD(:,sdftm<0),2));
        end
        if flag_halfwaverectify
            CSD(CSD > 0) = 0;
        end
        CSD = padarray(CSD,[1 0],NaN,'replicate');
    
        % Create Fig
         csdfig = figure('Position',[292 260 1042 550]);
         titleTextCSD = {'CSD',BRdatafile};
         sgtitle(titleTextCSD,'interpreter','none')
    
    
        % Subplot 1
        subplot(1,2,1)
        f_ShadedLinePlotbyDepth(CSD,corticaldepth,sdftm,[],1)
        set(gcf,'Position',[9.8000 49 510.4000 728.8000]); 
       
        % filterCSD
        CSDf = filterCSD(CSD);
        
        % Plot CSD
        subplot(1,2,2)
        switch sortdirection
            case 'ascending'
                y = chans;
                ydir = 'reverse';
            case 'descending'
                y = fliplr(chans);
                ydir = 'normal';
        end
        imagesc(sdftm,y,CSDf); colormap(flipud(jet));
        climit = max(abs(get(gca,'CLim'))*.8);
        set(gca,'CLim',[-climit climit],'Ydir',ydir,'Box','off','TickDir','out')
        hold on;
        plot([0 0], ylim,'k')
        c = colorbar;
        % caxis([-250 250])
        
        % Save figs
        saveFigNameFIG = strcat(OUTDIR,filesep,BRdatafile(23:end),'_CSD.fig');
        saveFigNamePNG = strcat(OUTDIR,filesep,BRdatafile(23:end),'_CSD.png');
        savefig(csdfig,saveFigNameFIG)
        saveas(csdfig,saveFigNamePNG)
        
    
    
    
    
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SECTION 3 -- PSD processing and plotting %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        % FFT
        Fs       = 1000; % Hz
        chanN    = size(SDF,1); % lfp is 11272360x24, SDF is ch x time x trial
        % loop through channels 
        for ch = 1:chanN
            %loop through trials
            for trialNum = 1:size(SDF,3)
                clear x n Spec
                lfp_holder        = SDF(ch,:,trialNum)';
            
                %%%%%%%%%%%%%%%%%
                % % % % % % % % % %  x =bandStopFiltLFP(lfp_holder); fix bandstop filter
                %%%%%%%%%%%%%%%%%
            
                n        = size(lfp_holder,2); % Number of data points
                % prep for psd 
                    nfft     = 512; 
                    window   = hanning(nfft); 
                    nwind    = length(window); 
                if n < nwind    % zero-pad x if it has length less than the window length
                    x(nwind)=0;  
                    n=nwind;
                end
                noverlap = 1;
                k        = fix((n-noverlap)/(nwind-noverlap));	% Number of windows
                index    = 1:nwind;
                % compute PSD
                    x        = lfp_holder;   
                    Spec     = zeros(nfft,1); 
                % loop through windows 
                    for j=1:k
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
                    power = nan(chanN,size(Spec,1),size(SDF,3));
                end
                power(ch,:,trialNum) = Spec; 
            end
        end
        % permute power, initial design was (spec x chanN)
        powerPermute = permute(power,[2 1 3]);
        
        % average power
        powerAvg = nanmean(powerPermute,3);
        
        
        % cheat at your band stop filter - remove 60Hz artifact manually
        idx60hz = find((freq_vector > 57.4 & freq_vector < 62.6 ));
        powerAvg(idx60hz,:) = 0;
        %remove top channel artifact?
            % warning('removing top channel artifact - necessary? BMC DEV!!')
            % powerAvg(:,1) = 0;
        
        % normalize power @ each frequency relative to power across contacts 
         power_norm = nan(size(powerAvg)); 
         for ch = 1:size(powerAvg,2)
             for f = 1:size(powerAvg,1)
                 power_norm(f,ch) = (powerAvg(f,ch) - mean(powerAvg(f,:)))./(mean(powerAvg(f,:))) * 100; % percent deviation from mean
             end
         end
        
    
    
        % Create Fig
            psdfig = figure('Position',[637 44.2000 749.6000 728.8000]);
            titleTextPSD = {'Stimulus Evoked PSD',BRdatafile};
            sgtitle(titleTextPSD,'interpreter','none')
    
        chans = 1:size(power_norm,2);
    
        subplot(1,2,1)
        set(gcf,'color','w'); 
        imagesc(freq_vector,chans,power_norm'); 
        colormap('hot'); xlim([0 100]); 
        xlabel('freq (Hz)'); ylabel('contact number'); 
        set(gca,'tickdir','out','ytick',chans); 
    
        
        % Get the Gamma x Beta cross
        % Beta is 12 - 20Hz (for our purposes)
        % Gamma is 30-59,61:100
        beta_index = (freq_vector > 12) & (freq_vector < 25);
        gamma_index = (freq_vector > 30) ;
        gamma_index(idx60hz) = false;
        
        clear j avgBeta avgGamma
        for j = chans
            avgBeta(j,1) = mean(power_norm(beta_index,j));
            avgGamma(j,1) = nanmean(power_norm(gamma_index,j));
        end
        
        
        % plot Gamma Beta Cross
        subplot(1,2,2)
        plot(avgBeta)
        hold on
        plot(fliplr(avgGamma))
        view([90 -90])
        set(gca,'xdir','reverse')
        legend('Beta','Gamma','Location','best')
        
        % Save figs
        saveFigNameFIG = strcat(OUTDIR,filesep,BRdatafile(23:end),'_PSD.fig');
        saveFigNamePNG = strcat(OUTDIR,filesep,BRdatafile(23:end),'_PSDprofile.png');
        savefig(psdfig,saveFigNameFIG)
        saveas(psdfig,saveFigNamePNG)
    
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SECTION 4 -- SAVE DATA                   %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Save data
        % avgBeta, avgGamma, power_norm, freq_vector, STIM, powerAvg. SDF800
        saveName = strcat(OUTDIR,filesep,BRdatafile(23:end),'_workspace.m');
        save(saveName,'-mat','-v7.3')
    
    
        %% End of filelist loop
        catch
            %nothing to do
        end  
    end
end