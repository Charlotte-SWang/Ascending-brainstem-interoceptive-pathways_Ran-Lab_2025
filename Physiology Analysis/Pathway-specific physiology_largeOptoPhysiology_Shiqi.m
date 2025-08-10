%% analyze physiology of the PVH and VLM stim experiment
% Ran 20201019
% Adapted by Shiqi 20250116
%clear;close all;clc;

[pathname_function] = uigetdir( 'Select the physiology function');
addpath(pathname_function)

%% Input

%%% Attention !!!%%%
SampleRate = 1; % default unit: ms --check 'isi' & 'isi_units'

% Read instructions for protocol extraction
[ProtocolName,ProtocolPath]=uigetfile('*.csv');
cd(ProtocolPath);
protocol_set = readmatrix(ProtocolName);

%%% Need to check
FixInfo = 6 % Channel is from 7:end

% channels
% Modified for Small Volume Injection 2025/01
ChannelName = ['ECG','Blood Pressure','Gastric Pressure','Breathing','EMG']
BChannel = 2;  
EKGChannel = 1;
GPChannel = 4;
BPChannel = 3;
laserChannel = 5;
%EMGChannel = 5


%% Analysis
ci = [1,2,3];
[fileDir] = uigetdir( 'Output');
cd(fileDir)
for iii = 3
    protocol_index=ci(iii);
    all_channel = find(protocol_set(ci(iii),FixInfo+1:end)>0)


    for jj = 1:(length(all_channel)-1)
        data_current = data_extract{protocol_index,all_channel(jj)};
        info_current = trialinfo_extract{protocol_index,all_channel(jj)};
        [unique_vals, ~, idx] = unique([info_current{:,1}],'stable');
        counts = accumarray(idx, 1);  % returns 5
        % Input parameter
        timeStart = 0
        timeEnd = protocol_set(protocol_index,4)*2.05;
        if iii == 3
            timeEnd =  protocol_set(protocol_index,4)*3;
        end
        
        sampleRate = 1000/SampleRate;
        % sampleRate = sampleRate / binN;
        if iii==2
            timePeriod = timeStart*sampleRate:timeEnd*sampleRate;
        else
            timePeriod = timeStart*sampleRate:timeEnd*sampleRate-1;
        end
        T = 1/sampleRate;
        L = length(timePeriod);             % Length of signal
        t = (0:L-1)*T;
        %nrep = 3


                         % laser onset and offset
        % laserData = data(timePeriod,laserChannel);
        % minLaserDistance = minLaserDistance * sampleRate;
        % laserData(laserData >= 4) = 4;
        % [plaseron,onsets] = findpeaks(laserData,timePeriod,'MinPeakDistance',minLaserDistance,'MinPeakHeight',2);
        % if length(onsets) ~= 1                  % only one trial
        %     laserDuration = (onsets(2) - onsets(1))/2;
        % else
        %     error('manual input of laser duration is required');
        laserDuration = protocol_set(protocol_index,4)*sampleRate;                     % Laser duration for 5 min stim
        
        % laserDuration = floor(laserDuration);
        % 
        % offsets = onsets + laserDuration;
        % figure;
        % plot(timePeriod,laserData,onsets,plaseron,'ro',offsets,plaseron,'ko');
        % % close all;
        % 
        % if length(onsets) > 10
        %     error('Check minLaserDistance, sampleRate, and timePeriod!');
        % end
        % 
        % minLength = min(offsets-onsets);        % make sure the time periods analyzed has equal length
        % lengthToAnalyze = floor(minLength*windowOfInterest); % length of the time periods to be analyzed
       
        % nTrial = length(min(onsets,offsets)); % number of trials to analyze
        BPwindowOfInterest = 1/6
        BPstimcalstart = 1/4
        BPlengthToAnalyze =  floor(laserDuration*BPwindowOfInterest);
        BPanalysisWindow = [];
        BPanalysisWindow(:,1) =[ (laserDuration*(1+BPstimcalstart)+1):1:(laserDuration*(1+BPstimcalstart)+BPlengthToAnalyze)];  % periods during stimulation to be analyzed
        BPanalysisWindow(:,2) =[ (laserDuration-BPlengthToAnalyze+1):1:laserDuration];       % periods before stimulation to be analyzed
   
        windowOfInterest = 1  %0.75 %use 1 for 5s ECG
        lengthToAnalyze =  floor(laserDuration*windowOfInterest);
        analysisWindow = [];
        analysisWindow(:,1) =[ (laserDuration*2-lengthToAnalyze+1):1:(laserDuration*2)];  % periods during stimulation to be analyzed
        analysisWindow(:,2) =[ (laserDuration-lengthToAnalyze+1): 1:laserDuration];       % periods before stimulation to be analyzed
        
        GPwindowOfInterest = 0.75 %0.75 %use 1 for 5s ECG
         GPlengthToAnalyze =  floor(laserDuration* GPwindowOfInterest);
         GPanalysisWindow = [];
         GPanalysisWindow(:,1) =[ (laserDuration*2- GPlengthToAnalyze+1):1:( laserDuration*2)];  % periods during stimulation to be analyzed
         GPanalysisWindow(:,2) =[ (laserDuration- GPlengthToAnalyze+1): 1: laserDuration];       % periods before stimulation to be analyzed
   
        BwindowOfInterest = 1
        BlengthToAnalyze =  floor(laserDuration*BwindowOfInterest);
        BanalysisWindow = [];
        BanalysisWindow(:,1) =[ (laserDuration+1):1:(laserDuration+BlengthToAnalyze)];  % periods during stimulation to be analyzed
        BanalysisWindow(:,2) =[ (laserDuration-BlengthToAnalyze+1): 1:laserDuration];       % periods before stimulation to be analyzed
        
        SwindowOfInterest = 1
        SlengthToAnalyze =  floor(laserDuration*SwindowOfInterest);
        SanalysisWindow = [];
        SanalysisWindow(:,1) =[ (laserDuration+1):1:(laserDuration+SlengthToAnalyze)];  % periods during stimulation to be analyzed
        SanalysisWindow(:,2) =[ (laserDuration-SlengthToAnalyze+1): 1:laserDuration];       % periods before stimulation to be analyzed
   
        % % %
 % % % % 
 %        %%%%%%%%%  Breathing-SIGH %%%%%%%%%%%%%
 %        if BChannel == all_channel(jj)
 % 
 %            %breathing parameters
 %            [FileSummary,FileSummaryPath]=uigetfile('*.xlsx');
 %            Sigh_parameter = readtable(strcat(FileSummaryPath,FileSummary),ReadRowNames=true);
 % 
 %            for k = 1:size(counts,1) %k = 1:size(info,1)
 %                if k<0 %
 %                    ifplot=1;
 %                else
 %                    ifplot=0;
 %                end
 % 
 %                for nn = 1:counts(k)
 %                    if k ==1
 %                        kk=nn;
 %                    else
 %                        kk = sum(counts(1:k-1)) +nn;
 %                    end
 % 
 %                    maxSighRate = 6 % estimated maximal breath rate in Hz;
 %                    ifInvert = Sigh_parameter{kk,1};      % recom: 0.18;
 %                    pSighThd =Sigh_parameter{kk,2};
 % 
 %                    breathData = detrend(data_current(:,kk));
 %                    for j = 1:2         % stim v. no stim
 %                        BreathingInWindow = breathData(SanalysisWindow(:,j));
 %                        [nSigh(kk,j)] = analyzeSigh(ifplot,protocol_index, k, nn,BreathingInWindow,SanalysisWindow(:,j),sampleRate,maxSighRate,pSighThd,ifInvert);     
 %                    end   
 % 
 %                    changeNSigh(kk) = nSigh(kk,1) - nSigh(kk,2);
 %                end
 %                for j =1:2
 %                    nSigh_mouse(k,j) = mean(nSigh(sum(counts(1:k-1))+1:sum(counts(1:k)),j));
 %                end
 % 
 %                  changeNSigh_mouse(k) = mean(changeNSigh(sum(counts(1:k-1))+1:sum(counts(1:k))));
 % 
 %            end
 % 
 %           % %%% write data
 %            clear writeBR; clear writeBR_mouse;
 %            writeBR = {'Order','Mouse #','PVH/VLM','C/E','Trial #',...
 %                'nSigh Post Stim','nSigh Pre Stim','nSigh Change #'};
 %            for i = 1:sum(counts)
 %                writeBR = [writeBR;{info_current{i,1},info_current{i,2},info_current{i,3},info_current{i,4},info_current{i,5},...
 %                   nSigh(i,1),nSigh(i,2),changeNSigh(i)}];
 %            end
 %            writeBR = [writeBR;{[],[],[],timeStart,timeEnd,sampleRate,'pSighThd',pSighThd}];
 % 
 %            writeBR_mouse = {'Order','Mouse #','PVH/VLM','C/E',...
 %                'nSigh Post Stim','nSigh Pre Stim','nSigh Change #'};
 %            for i = 1:size(counts,1)
 %                writeBR_mouse = [writeBR_mouse;{info_current{sum(counts(1:i)),1},info_current{sum(counts(1:i)),2},info_current{sum(counts(1:i)),3},info_current{sum(counts(1:i)),4},...
 %                    nSigh_mouse(i,1),nSigh_mouse(i,2),changeNSigh_mouse(i)}];
 %            end
 %            writeBR_mouse = [writeBR_mouse;{[],[],timeStart,timeEnd,sampleRate,'pSighThd',pSighThd}];
 % 
 %            writecell(writeBR,['Protocol',num2str(protocol_index),'_Sigh_trial.csv']);
 %            writecell(writeBR_mouse,['Protocol',num2str(protocol_index),'_Sigh_mouse.csv']);
 %        end
 %        % 

        %%%%%%%%%  Breathing %%%%%%%%%%%%%
        if BChannel == all_channel(jj)

            %breathing parameters
            [FileSummary,FileSummaryPath]=uigetfile('*.xlsx');
            Breathing_parameter = readtable(strcat(FileSummaryPath,FileSummary),ReadRowNames=true);
            BreathRate = {};
            amp_stage = {};
            amp_stage_tv = {};
            BR_zscore = {};
            amp_window = {};
            amp_stage_raw = {};
            BR_CI={};
            Delta_amp = {}
            CI_amp = {};
            Delta_amp_tv={};
            CI_amp_tv={};
            single_ci_amp_tv={};
            single_z_amp_tv={};
            single_z_amp_tv_mouse={};
            single_z_amp_tv_trial={};
            TV_auc_ci={};


            %for k = 34% 
            for k = 1:size(counts,1) %k = 1:size(info,1)
                if k<0 %
                %if k >=30 & k <=34 %<0 %44 
                    ifplot=1;
                else
                    ifplot=0;
                end

                for nn = 1:counts(k)
                    if k ==1
                        kk=nn;
                    else
                        kk = sum(counts(1:k-1)) +nn;
                    end

                    maxBreathRate = Breathing_parameter{kk,1};  
                    %maxBreathRate = 6 % estimated maximal breath rate in Hz;
                    %pbreaththresholdlow = 0.20;      % recom: 0.18;
                    %pbreaththreshold = 0.8;         % recom: 0.5; (the lower the number, the higher the thd)
                    pbreaththresholdlow = Breathing_parameter{kk,3};      % recom: 0.18;
                    pbreaththreshold =Breathing_parameter{kk,2};
                    pSighThd = Breathing_parameter{kk,4};% recom: 0.5; (the lower the number, the higher the thd)

                    breathData = detrend(data_current(:,kk));
                    for j = 1:2         % stim v. no stim
                        BreathingInWindow = breathData(BanalysisWindow(:,j));
                        [breathRate(kk,j),mVolume(kk,j),TVolume(kk,j),TVolume_auc(kk,j),TVolume_auc_nogasp(kk,j),mVolume_nogasp(kk,j),nSigh(kk,j), abspbthreshold,abspbthresholdlow,breathEvents] = analyzeBreath(ifplot,protocol_index, k, nn,BreathingInWindow,BanalysisWindow(:,j),sampleRate,maxBreathRate,pbreaththreshold,pbreaththresholdlow,pSighThd);     
                    end   
                    changeBreathRate(kk) = (breathRate(kk,1) - breathRate(kk,2))/(breathRate(kk,1) +breathRate(kk,2));
                    changeTVolume(kk) = (TVolume(kk,1) - TVolume(kk,2))/(TVolume(kk,1)+TVolume(kk,2));   
                    %%changeNSigh(i) = (nSigh(i,1) - nSigh(i,2))/nSigh(i,2) * 100;
                  
                    changeTVolume_auc(kk) = (TVolume_auc(kk,1) - TVolume_auc(kk,2))/(TVolume_auc(kk,1)+TVolume_auc(kk,2));   
                    changemVolume(kk) =(mVolume(kk,1) - mVolume(kk,2))/(mVolume(kk,1)+mVolume(kk,2));    
                    changeTVolume_auc_nogasp(kk) = (TVolume_auc_nogasp(kk,1) - TVolume_auc_nogasp(kk,2))/(TVolume_auc_nogasp(kk,1)+TVolume_auc_nogasp(kk,2));   
                    changemVolume_nogasp(kk) =(mVolume_nogasp(kk,1) - mVolume_nogasp(kk,2))/(mVolume_nogasp(kk,1)+mVolume_nogasp(kk,2));    
                    
                    changeNSigh(kk) = nSigh(kk,1) - nSigh(kk,2);

                    nnbin = 5; % Change!!! number of bins in each stage
                    bin = laserDuration/nnbin; % bin
                    tt = timeEnd*sampleRate;
                    bintime = bin/sampleRate;
                    baseline = [1:nnbin];
                    baseline = [nnbin-4:nnbin];

                    baseline_ms=BanalysisWindow(:,1);


                    [TV_auc_ci{kk},TV_auc_ci_nogasp{kk},single_z_amp_tv{kk},single_ci_amp_tv{kk},BR_CI{kk}, Delta_amp{kk}, CI_amp{kk}, Delta_amp_tv{kk}, CI_amp_tv{kk},BreathRate{kk},BR_zscore{kk},amp_stage{kk},amp_stage_tv{kk},amp_window{kk},amp_stage_raw{kk}] = analyzeBreath_SW(BlengthToAnalyze,protocol_index, k, nn,tt,bin,baseline,baseline_ms,laserDuration,breathData,timePeriod,sampleRate,maxBreathRate,pbreaththreshold,pbreaththresholdlow,pSighThd, abspbthreshold,abspbthresholdlow);
                    changeMinuteVolume(kk) = (sum(amp_window{kk}{:,1}) - sum(amp_window{kk}{:,2}))/(sum(amp_window{kk}{:,1})+sum(amp_window{kk}{:,2}));
                    single_z_amp_tv_trial{kk}= mean(single_z_amp_tv{kk});
                
                end
                for j =1:2
                    breathRate_mouse(k,j) = mean(breathRate(sum(counts(1:k-1))+1:sum(counts(1:k)),j));
                    TVolume_mouse(k,j) = mean(TVolume(sum(counts(1:k-1))+1:sum(counts(1:k)),j));
                    
                    TVolume_auc_mouse(k,j) = mean(TVolume_auc(sum(counts(1:k-1))+1:sum(counts(1:k)),j));
                    mVolume_mouse(k,j) = mean(mVolume(sum(counts(1:k-1))+1:sum(counts(1:k)),j));
                    TVolume_auc_nogasp_mouse(k,j) = mean(TVolume_auc_nogasp(sum(counts(1:k-1))+1:sum(counts(1:k)),j));
                    mVolume_nogasp_mouse(k,j) = mean(mVolume_nogasp(sum(counts(1:k-1))+1:sum(counts(1:k)),j));
                    
                    nSigh_mouse(k,j) = mean(nSigh(sum(counts(1:k-1))+1:sum(counts(1:k)),j));
                end

                    single_z_amp_tv_mouse{k}= mean(cat(1,single_z_amp_tv{sum(counts(1:k-1))+1:sum(counts(1:k))}));
                    changeBreathRate_mouse(k) = mean(changeBreathRate(sum(counts(1:k-1))+1:sum(counts(1:k))));
                    changeTVolume_mouse(k) = mean(changeTVolume(sum(counts(1:k-1))+1:sum(counts(1:k))));
                   
                    changeTVolume_auc_mouse(k) = mean(changeTVolume_auc(sum(counts(1:k-1))+1:sum(counts(1:k))));
                    changemVolume_mouse(k) = mean(changemVolume(sum(counts(1:k-1))+1:sum(counts(1:k))));
                    changeTVolume_auc_nogasp_mouse(k) = mean(changeTVolume_auc_nogasp(sum(counts(1:k-1))+1:sum(counts(1:k))));
                    changemVolume_nogasp_mouse(k) = mean(changemVolume_nogasp(sum(counts(1:k-1))+1:sum(counts(1:k))));
                    
                    changeNSigh_mouse(k) = mean(changeNSigh(sum(counts(1:k-1))+1:sum(counts(1:k))));
                    changeMinuteVolume_mouse(k) = mean(changeMinuteVolume(sum(counts(1:k-1))+1:sum(counts(1:k))));
              end

            %%%%%%write data
            clear writeBR; clear writeBR_mouse;
            writeBR = {'Order','Mouse #',...
                'PVH/VLM','C/E','Trial #',...
                'Breath Rate Post Stim','Breath Rate Pre Stim','Breath Rate Change %',...
                'TV Post Stim','TV Pre Stim','TV Change %',...
                'TV AUC Post Stim','TV AUC Pre Stim','TV AUC Change %',...
                'MV Post Stim','MV Pre Stim','MV Change %',...
                'TV AUC nogasp Post Stim','TV AUC nogasp Pre Stim','TV AUC nogasp Change %',...
                'MV nogasp Post Stim','MV nogasp Pre Stim','MV nogasp Change %',...
                'nSigh Post Stim','nSigh Pre Stim','nSigh Change #',...
                'Minute Volume Change %'};
            for i = 1:sum(counts)
                writeBR = [writeBR;{info_current{i,1},info_current{i,2},...
                    info_current{i,3},info_current{i,4},info_current{i,5},...
                    breathRate(i,1),breathRate(i,2),changeBreathRate(i),...
                    TVolume(i,1),TVolume(i,2),changeTVolume(i),...
                    TVolume_auc(i,1),TVolume_auc(i,2),changeTVolume_auc(i),...
                    mVolume(i,1),mVolume(i,2),changemVolume(i),...
                    TVolume_auc_nogasp(i,1),TVolume_auc_nogasp(i,2),changeTVolume_auc_nogasp(i),...
                    mVolume_nogasp(i,1),mVolume_nogasp(i,2),changemVolume_nogasp(i),...
                    nSigh(i,1),nSigh(i,2),changeNSigh(i),changeMinuteVolume(i)}];
            end
            writeBR = [writeBR;{[],[],...
                [],[],[],...
                timeStart,timeEnd,sampleRate,...
                'pSighThd',pSighThd,[],...
                [],[],[],...
                [],[],[],...
                [],[],[],...
                [],[],[],...
                [],[],[],[]}];

            writeBR_mouse = {'Order','Mouse #','PVH/VLM','C/E',...
                'Breath Rate Post Stim','Breath Rate Pre Stim','Breath Rate Change %',...
                'TV Post Stim','TV Pre Stim','TV Change %',...
                 'TV AUC Post Stim','TV AUC Pre Stim','TV AUC Change %',...
                'MV Post Stim','MV Pre Stim','MV Change %',...
                'TV AUC nogasp Post Stim','TV AUC nogasp Pre Stim','TV AUC nogasp Change %',...
                'MV nogasp Post Stim','MV nogasp Pre Stim','MV nogasp Change %',...
                'nSigh Post Stim','nSigh Pre Stim','nSigh Change #',...
                'Minute Volume Change %'};
            for i = 1:size(counts,1)
                writeBR_mouse = [writeBR_mouse;{info_current{sum(counts(1:i)),1},info_current{sum(counts(1:i)),2},...
                    info_current{sum(counts(1:i)),3},info_current{sum(counts(1:i)),4},...
                    breathRate_mouse(i,1),breathRate_mouse(i,2),changeBreathRate_mouse(i),...
                    TVolume_mouse(i,1),TVolume_mouse(i,2),changeTVolume_mouse(i),...
                    TVolume_auc_mouse(i,1),TVolume_auc_mouse(i,2),changeTVolume_auc_mouse(i),...
                    mVolume_mouse(i,1),mVolume_mouse(i,2),changemVolume_mouse(i),...
                    TVolume_auc_nogasp_mouse(i,1),TVolume_auc_nogasp_mouse(i,2),changeTVolume_auc_nogasp_mouse(i),...
                    mVolume_nogasp_mouse(i,1),mVolume_nogasp_mouse(i,2),changemVolume_nogasp_mouse(i),...
                    nSigh_mouse(i,1),nSigh_mouse(i,2),changeNSigh_mouse(i),changeMinuteVolume_mouse(i)}];
            end
            writeBR_mouse = [writeBR_mouse;{[],[],...
                [],[],...
                timeStart,timeEnd,sampleRate,...
                'pSighThd',pSighThd,[],...
                [],[],[],...
                [],[],[],...
                [],[],[],...
                [],[],[],...
                [],[],[],[]}];

            writecell(writeBR,['250806_Protocol',num2str(protocol_index),'_Breathing_trial.csv']);
            writecell(writeBR_mouse,['250806_Protocol',num2str(protocol_index),'_Breathing_mouse.csv']);

            matname =['250806_Protocol',num2str(protocol_index),'_Breathing_plot_data.mat']
            save (matname,'trialinfo_extract','TV_auc_ci','TV_auc_ci_nogasp','single_z_amp_tv_mouse','single_z_amp_tv_trial','single_z_amp_tv','single_ci_amp_tv','BR_zscore','BR_CI','protocol_index', 'BreathRate', 'amp_stage','Delta_amp','CI_amp','amp_stage_tv','Delta_amp_tv','CI_amp_tv','amp_stage_raw','bin')
            bintime = bin/sampleRate
        end

        % %%%%%%%%%  EKG %%%%%%%%%%%%%
        % if EKGChannel == all_channel(jj)
        % 
        %     maxHeartRate = 15;      % estimated maximal breath rate in Hz;
        %     EKGThd = -10;
        %     EKGThdLow = 10;
        % 
        %     HeartRate = {};
        %     amp_stage = {};
        %     HR_zscore = {};
        %     for k = 1: size(counts,1)
        % 
        %         %%%%%% Tune parameters for 5s heart rate--Start
        %         if iii ==3 
        %             if k== 29
        %                 maxHeartRate = 8;
        %             elseif k == 22
        %                 maxHeartRate = 10;
        %             elseif k==34
        %                 maxHeartRate = 13;
        %             else
        %                 maxHeartRate = 12;
        %             end
        %             if k== 34 |k==29
        %                 EKGThd = -0.2;
        %                 EKGThdLow = -0.9;
        %             elseif k ==22
        %                 EKGThd = 0.07;
        %                 EKGThdLow = -0.4;
        %             elseif k== 10|k== 11| k== 12| k==13| k==7|k==31|k==32
        %                 EKGThd = -10;
        %                 EKGThdLow = 10;
        %             else
        %                 EKGThd = 0.15;
        %                 EKGThdLow =-0.15;
        %             end
        %         end
        %         %%%%%% Tune parameters for 5s heart rate--End
        % 
        %         for nn = 1:counts(k)
        %             if k ==1
        %                 kk=nn;
        %             else
        %                 kk = sum(counts(1:k-1)) +nn;
        %             end
        %             % EKG parameters
        %             % Constant value for now
        % 
        %             EKGData = data_current(:,kk);
        %             for j = 1:2         % stim v. no stim
        %                 EKGInWindow = EKGData(analysisWindow(:,j));
        %                 [heartRate(kk,j),EKGAmplitude(kk,j)] = analyzeEKG(protocol_index, k,EKGInWindow,analysisWindow(:,j),sampleRate,maxHeartRate,EKGThd,EKGThdLow);
        %             end
        %             changeHeartRate(kk) = (heartRate(kk,1) - heartRate(kk,2))/(heartRate(kk,1) + heartRate(kk,2));
        %             changeEKGAmplitude(kk) = (EKGAmplitude(kk,1) - EKGAmplitude(kk,2))/(EKGAmplitude(kk,1) + EKGAmplitude(kk,2));
        % 
        %             nbin = 60; % Change!!! number of bins in each stage
        %             bin = laserDuration/nbin; % ms
        %             tt = timeEnd*sampleRate; % ms
        %             bintime = bin/sampleRate;
        %             baseline = [1:nbin];
        %             [HeartRate{kk},HR_zscore{kk},amp_stage{kk},HR_CI{kk}] = analyzeEKG_SW(protocol_index, k, nn,tt,bin,baseline,laserDuration,EKGData,timePeriod,sampleRate,maxHeartRate,EKGThd,EKGThdLow);
        % 
        %         end
        %         for j =1:2
        %             heartRate_mouse(k,j) = mean(heartRate(sum(counts(1:k-1))+1:sum(counts(1:k)),j));
        %             EKGAmplitude_mouse(k,j) = mean(EKGAmplitude(sum(counts(1:k-1))+1:sum(counts(1:k)),j));
        %         end
        %         changeHeartRate_mouse(k) = mean(changeHeartRate(sum(counts(1:k-1))+1:sum(counts(1:k))));
        %         changeEKGAmplitude_mouse(k) = mean(changeEKGAmplitude(sum(counts(1:k-1))+1:sum(counts(1:k))));
        %     end
        % 
        %     % write data
        %     clear writeEKG; clear writeEKG_mouse;
        %     writeEKG = {'Order','Mouse #','PVH/VLM','C/E','Trial #','Heart Rate Post Stim','Heart Rate Pre Stim','Heart Rate Change %','EKG Amp Post Stim','EKG Amp Pre Stim','EKG Amp Change %'};
        %     for i = 1:sum(counts)
        %         writeEKG = [writeEKG;{info_current{i,1},info_current{i,2},info_current{i,3},info_current{i,4},info_current{i,5},heartRate(i,1),heartRate(i,2),changeHeartRate(i),EKGAmplitude(i,1),EKGAmplitude(i,2),changeEKGAmplitude(i)}];
        %     end
        %     writeEKG = [writeEKG;{[],[],[],[],[],timeStart,timeEnd,sampleRate,'EKGThd',EKGThd,[]}];
        % 
        %     writeEKG_mouse = {'Order','Mouse #','PVH/VLM','C/E','Heart Rate Post Stim','Heart Rate Pre Stim','Heart Rate Change %','EKG Amp Post Stim','EKG Amp Pre Stim','EKG Amp Change %'};
        %     for i = 1:size(counts,1)
        %         writeEKG_mouse = [writeEKG_mouse;{info_current{sum(counts(1:i)),1},info_current{sum(counts(1:i)),2},info_current{sum(counts(1:i)),3},info_current{sum(counts(1:i)),4},heartRate_mouse(i,1),heartRate_mouse(i,2),changeHeartRate_mouse(i),EKGAmplitude_mouse(i,1),EKGAmplitude_mouse(i,2),changeEKGAmplitude_mouse(i)}];
        %     end
        %     writeEKG_mouse = [writeEKG_mouse;{[],[],[],[],timeStart,timeEnd,sampleRate,'EKGThd',EKGThd,[]}];
        % 
        %     %[fileDir] = uigetdir( 'Output');
        %     %cd(fileDir);
        %     writecell(writeEKG,['Protocol',num2str(protocol_index),'_HeartRate_trial.csv']);
        %     writecell(writeEKG_mouse,['Protocol',num2str(protocol_index),'_HeartRate_mouse.csv']);
        % 
        %     matname =['Protocol',num2str(protocol_index),'_HeartRate_plot_data.mat']
        %     save (matname,'HR_CI','HR_zscore','protocol_index', 'HeartRate', 'amp_stage')
        % 
        % end


        % %%%%%%%%  Blood Pressure %%%%%%%%%%%%%
        % if BPChannel == all_channel(jj)
        %     EstMaxBP = 11; %Hz
        % 
        %     MinPeakH_max = -10;
        %     MinPeakH_min = 10;
        % 
        %     SBP = {};
        %     SBP_zscore = {};
        % 
        %     BP = {};
        %     BP_zscore = {};
        % 
        %     DBP = {};
        %     DBP_zscore = {};
        % 
        %     MAP = {};
        %     MAP_zscore = {};
        %     CI_BP ={}
        %     Delta_BP ={}
        % 
        %     for k = 1:size(counts,1) %k = 1:size(info,1)
        %         if k == 9
        %             EstMaxBP = 20;
        %         else 
        %             EstMaxBP = 11;
        %         end
        %         if k<0 %k >=34 & k <=35 %<0 %44 
        %             ifplot=1;
        %         else
        %             ifplot=0;
        %         end
        % 
        %         for nn = 1:counts(k)
        %             if k ==1
        %                 kk=nn
        %             else
        %                 kk = sum(counts(1:k-1)) +nn
        %             end
        % 
        %             BPData = data_current(:,kk);
        %             figure
        %             plot(BPData)
        %             hold on
        %             plot(bandstop(BPData,[0.2 3],sampleRate))
        %             BPData = bandstop(BPData,[0.2 3],sampleRate);
        %             for j = 1:2         % stim v. no stim
        %                 BPInWindow = BPData(BPanalysisWindow(:,j));
        %                 [BP_windowmean(kk,j),SBP_windowmean(kk,j),DBP_windowmean(kk,j),PulsePressure_windowmean(kk,j)] = analyzeBP(protocol_index, k, nn,BPInWindow,BPanalysisWindow(:,j),sampleRate,EstMaxBP,MinPeakH_max,MinPeakH_min);
        %             end   
        %             changeSBP(kk) = (SBP_windowmean(kk,1) - SBP_windowmean(kk,2))/(SBP_windowmean(kk,1) + SBP_windowmean(kk,2));
        %             changeDBP(kk) = (DBP_windowmean(kk,1) - DBP_windowmean(kk,2))/(DBP_windowmean(kk,1) + DBP_windowmean(kk,2));
        %             changePP(kk) = (PulsePressure_windowmean(kk,1) - PulsePressure_windowmean(kk,2))/(PulsePressure_windowmean(kk,1) + PulsePressure_windowmean(kk,2));
        %             deltaBP(kk)=(SBP_windowmean(kk,1)+DBP_windowmean(kk,1))/2-(SBP_windowmean(kk,2)+DBP_windowmean(kk,2))/2;
        % 
        %             nbin = 60; % Change!!! number of bins in stimulation period
        %             bin = laserDuration/nbin; % ms
        %             tt = timeEnd*1000; % ms
        %             baseline = [nbin-10:nbin];
        % 
        %             [Delta_BP{kk}, CI_BP{kk},SBP{kk},SBP_zscore{kk},DBP{kk},DBP_zscore{kk},MAP{kk},MAP_zscore{kk}] = analyzeBP_SW(protocol_index, k, nn,tt,bin,baseline,laserDuration,BPData,timePeriod,sampleRate,EstMaxBP,MinPeakH_max,MinPeakH_min);
        % 
        %         end
        %         for j =1:2
        %             SBP_mouse(k,j) = mean(SBP_windowmean(sum(counts(1:k-1))+1:sum(counts(1:k)),j));
        %             DBP_mouse(k,j) = mean(DBP_windowmean(sum(counts(1:k-1))+1:sum(counts(1:k)),j));
        %             PulsePressure_mouse(k,j) = mean(PulsePressure_windowmean(sum(counts(1:k-1))+1:sum(counts(1:k)),j));
        %         end
        %             changeSBP_mouse(k) = mean(changeSBP(sum(counts(1:k-1))+1:sum(counts(1:k))));
        %             changeDBP_mouse(k) = mean(changeDBP(sum(counts(1:k-1))+1:sum(counts(1:k))));
        %             changePP_mouse(k) = mean(changePP(sum(counts(1:k-1))+1:sum(counts(1:k))));
        %             deltaBP_mouse(k) = mean(deltaBP(sum(counts(1:k-1))+1:sum(counts(1:k))));
        %     end
        % 
        %     %%write data
        %     clear writeBP; clear writeBP_mouse;
        %     writeBP = {'Order','Mouse #','PVH/VLM','C/E','Trial #','deltaBP','SBP Post Stim','SBP Pre Stim','SBP Change %','DBP Post Stim','DBP Pre Stim','DBP Change %','PulsePresure Post Stim','PulsePresure Pre Stim','PulsePresure Change %'};
        %     for i = 1:sum(counts)
        %         writeBP = [writeBP;{info_current{i,1},info_current{i,2},info_current{i,3},info_current{i,4},info_current{i,5},deltaBP,SBP_windowmean(i,1),SBP_windowmean(i,2),changeSBP(i),DBP_windowmean(i,1),DBP_windowmean(i,2),changeDBP(i),PulsePressure_windowmean(i,1),PulsePressure_windowmean(i,2),changePP(i)}];
        %     end
        %     writeBP = [writeBP;{[],[],[],[],[],timeStart,timeEnd,sampleRate,[],[],[],[],[],[],[]}];
        % 
        %     writeBP_mouse = {'Order','Mouse #','PVH/VLM','C/E','deltaBP','SBP Post Stim','SBP Pre Stim','SBP Change %','DBP Post Stim','DBP Pre Stim','DBP Change %','PulsePresure Post Stim','PulsePresure Pre Stim','PulsePresure Change %'};
        %     for i = 1:size(counts,1)
        %         writeBP_mouse = [writeBP_mouse;{info_current{sum(counts(1:i)),1},info_current{sum(counts(1:i)),2},info_current{sum(counts(1:i)),3},info_current{sum(counts(1:i)),4}, deltaBP_mouse(i),SBP_mouse(i,1),SBP_mouse(i,2),changeSBP_mouse(i),DBP_mouse(i,1),DBP_mouse(i,2),changeDBP_mouse(i),PulsePressure_mouse(i,1),PulsePressure_mouse(i,2),changePP_mouse(i)}];
        %     end
        %     writeBP_mouse = [writeBP_mouse;{[],[],[],[],timeStart,timeEnd,sampleRate,[],[],[],[],[],[],[]}];
        %     [fileDir] = uigetdir( 'Output');
        %     cd(fileDir);
        %     writecell(writeBP,['Protocol',num2str(protocol_index),'_BloodPressure_trial.csv']);
        %     writecell(writeBP_mouse,['Protocol',num2str(protocol_index),'_BloodPressure_mouse.csv']);
        % 
        %     matname =['Protocol',num2str(protocol_index),'_BloodPressure_delta_plot_data.mat'];
        %     save (matname,'Delta_BP','CI_BP','SBP','SBP_zscore','DBP','DBP_zscore','MAP','MAP_zscore','protocol_index');
        %     bintime = bin/sampleRate;
        % end


        % %%%%%%%%%%  GP %%%%%%%%%%%%%
        % if GPChannel == all_channel(jj)
        %     % GP parameters
        %     lowpassCutoffGP = 1.5;
        %     arbitraryBcgdPeriod = 5;
        %     nDetrend = 3;
        %     conflvl = 0.9;
        %     GP={};
        %     GP_Phasic ={};
        %     GP_Tonic = {};
        %     clear detrendGPdata trendGP
        %     for k = 1: size(counts,1)
        %         for nn = 1:counts(k)
        %             if k ==1
        %                 kk=nn
        %             else
        %                 kk = sum(counts(1:k-1)) +nn
        %             end
        % 
        %             samplepts = linspace(0,1, GPlengthToAnalyze);
        %             GPData = data_current(:,kk);
        %             GP{k,nn}= GPData;
        %             if isempty(GPData)
        %                 comtinue
        %             end
        %             samplepts_whole = linspace(0,1,length(GP{k,nn}));
        %             [polorder,GP_Phasic{k,nn},GP_Tonic{k,nn},] = RobustDetrend(GP{k,nn},nDetrend,conflvl,samplepts_whole);
        % 
        % 
        %             for j = 1:2         % stim v. no stim
        %                 GPInWindow = GPData( GPanalysisWindow(:,j));  
        %                 [polorder,detrendGPdata(:,j),trendGP(:,j)] = RobustDetrend(GPInWindow,nDetrend,conflvl,samplepts);
        %                 detrendGPdata(:,j) = smooth(detrendGPdata(:,j),sampleRate);
        % 
        %                 %figure;plot(squeeze(detrendGPdata(:,j)));hold on;plot(squeeze(trendGP(:,j)));hold on; plot(GPInWindow)
        %                 %title(sprintf('Animal #%d Trial #%d, Stim #%d', k,nn,j));
        % 
        %                 phasic(kk,j) = sum(detrendGPdata(:,j))/(300*sampleRate); %% Unit: cmH2O x min
        %                 tonic(kk,j) = mean(trendGP(:,j));
        %                 GPInWindow = smooth(GPInWindow);
        %                 % total(i,j) = sum(GPInWindow-min(GPInWindow)); 
        %                 total(kk,j) = sum(GPInWindow)/(5*sampleRate); %% Unit: cmH2O x min
        %                 total_mean(kk,j)=mean(GPInWindow);
        %             end  
        %             minGP = min(trendGP(:,:),[],'all')           
        %             for j = 1:2
        %                 GPInWindow = GPData( GPanalysisWindow(:,j)); 
        %                 total_sbs(kk,j) = sum(GPInWindow-minGP)
        %             end
        % 
        % 
        %             changePhasic(kk) = (phasic(kk,1) - phasic(kk,2))/(phasic(kk,1) + phasic(kk,2));
        %             changeTonic(kk) = (tonic(kk,1) - tonic(kk,2))/(tonic(kk,1) + tonic(kk,2));
        %             changeTotal(kk) = (total(kk,1) - total(kk,2))/(total(kk,1) + total(kk,2));
        %             changeTotal_mean(kk) = (total_mean(kk,1) - total_mean(kk,2))/ (total_mean(kk,1) + total_mean(kk,2));
        %             changeTotal_sbs(kk) = (total_sbs(kk,1) - total_sbs(kk,2))/(total_sbs(kk,1) + total_sbs(kk,2));
        %             deltaTotal(kk)= (total(kk,1) - total(kk,2));
        %         end
        % 
        %         for j =1:2
        %             phasic_mouse(k,j) = mean(phasic(sum(counts(1:k-1))+1:sum(counts(1:k)),j));
        %             tonic_mouse(k,j) = mean(tonic(sum(counts(1:k-1))+1:sum(counts(1:k)),j));
        %             total_mouse(k,j) = mean(total(sum(counts(1:k-1))+1:sum(counts(1:k)),j));
        %             total_mean_mouse(k,j) = mean(total_mean(sum(counts(1:k-1))+1:sum(counts(1:k)),j));
        %             total_sbs_mouse(k,j) =  mean(total_sbs(sum(counts(1:k-1))+1:sum(counts(1:k)),j));
        %         end
        % 
        %         changePhasic_mouse(k) = mean(changePhasic(sum(counts(1:k-1))+1:sum(counts(1:k))));
        %         changeTonic_mouse(k) = mean(changeTonic(sum(counts(1:k-1))+1:sum(counts(1:k))));
        %         changeTotal_mouse(k) = mean(changeTotal(sum(counts(1:k-1))+1:sum(counts(1:k))));
        %         changeTotal_mean_mouse(k) = mean(changeTotal_mean(sum(counts(1:k-1))+1:sum(counts(1:k))));
        %         changeTotal_sbs_mouse(k) = mean(changeTotal_sbs(sum(counts(1:k-1))+1:sum(counts(1:k))));
        %         deltaTotal_mouse(k)= mean(deltaTotal(sum(counts(1:k-1))+1:sum(counts(1:k))));
        %     end
        %     % write data
        %     clear writeGP;clear writeGP_mouse;
        %     writeGP = {'Order','Mouse #','PVH/VLM','C/E','Trial #',...
        %         'Phsic Post Stim','Phasic Pre Stim','Phsic Change Index',...
        %         'Tonic Post Stim','Tonic Pre Stim','Tonic Change Index',...
        %         'Total Post Stim','Total Pre Stim','Total Change Index',...
        %         'Total Subtract Post Stim','Total Subtract Pre Stim', 'Total Subtract Change Index',...
        %         'Delta Total'};
        %     for i = 1:sum(counts)
        %         writeGP = [writeGP;{info_current{i,1},info_current{i,2},...
        %             info_current{i,3},info_current{i,4},info_current{i,5},...
        %             phasic(i,1),phasic(i,2),changePhasic(i),...
        %             tonic(i,1),tonic(i,2),changeTonic(i),...
        %             total(i,1),total(i,2),changeTotal(i), total_sbs(i,1),total_sbs(i,2),changeTotal_sbs(i), deltaTotal(i)}];
        %     end
        %     writeGP = [writeGP;{[],[],[],[],[],timeStart,timeEnd,sampleRate,[],[],[],[],[],[],[],[],[],[]}];
        % 
        %     writeGP_mouse = {'Order','Mouse #','PVH/VLM','C/E',...
        %         'Phsic Post Stim','Phasic Pre Stim','Phsic Change Index',...
        %         'Tonic Post Stim','Tonic Pre Stim','Tonic Change Index',...
        %         'Total Post Stim','Total Pre Stim','Total Change Index',...
        %         'Total Subtract Post Stim','Total Subtract Pre Stim', 'Total Subtract Change Index',...
        %         'Delta Total'};
        %     for i = 1:size(counts,1)
        %         writeGP_mouse = [writeGP_mouse;{info_current{sum(counts(1:i)),1},info_current{sum(counts(1:i)),2},...
        %             info_current{sum(counts(1:i)),3},info_current{sum(counts(1:i)),4},...
        %             phasic_mouse(i,1),phasic_mouse(i,2),changePhasic_mouse(i),...
        %             tonic_mouse(i,1),tonic_mouse(i,2),changeTonic_mouse(i),...
        %             total_mouse(i,1),total_mouse(i,2),changeTotal_mouse(i), ...
        %             total_sbs_mouse(i,1),total_sbs_mouse(i,2),changeTotal_sbs_mouse(i), ...
        %             deltaTotal_mouse(i)}];
        %     end
        %     writeGP_mouse = [writeGP_mouse;{[],[],[],[],timeStart,timeEnd,sampleRate,[],[],[],[],[],[],[],[],[],[]}];
        %     % [fileDir] = uigetdir( 'Output');
        %     % cd(fileDir);
        %     writecell(writeGP,['Protocol',num2str(protocol_index),'_GastricPressure_trial.csv']);
        %     writecell(writeGP_mouse,['Protocol',num2str(protocol_index),'_GastricPressure_mouse.csv']);
        % 
        %     matname =['Protocol',num2str(protocol_index),'_GastricPressure_plot_data.mat'];
        %     save (matname,'GP','GP_Phasic','GP_Tonic','protocol_index');
        %     %bintime = bin/sampleRate;
        % end

    end
end

          
% %%%%%%%%%% Close hidden figures 
% invisFigs = findall(0, 'Type', 'figure', 'Visible', 'off');
% close(invisFigs);
   
%% %%%%%%%%%% Close all figures 
allFigs = findall(0, 'Type', 'figure');
allFigs(allFigs == gcf) = [];
close(allFigs);