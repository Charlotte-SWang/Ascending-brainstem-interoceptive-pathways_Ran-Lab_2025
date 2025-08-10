%% analyze physiology of the PVH and VLM stim experiment
% Ran 20201019
% Adapted by Shiqi 20250116
%clear;close all;clc;

[pathname_function] = uigetdir( 'Select the physiology function');
addpath(pathname_function)

%% Input

%%% Attention !!!%%%
SampleRate = 0.2; % default unit: ms --check 'isi' & 'isi_units'

% Read instructions for protocol extraction
[ProtocolName,ProtocolPath]=uigetfile('*.csv');
cd(ProtocolPath);
protocol_set = readmatrix(ProtocolName);

%%% Need to check
FixInfo = 6 % Channel is from 7:end

% channels
% Modified for Small Volume Injection 2025/01
ChannelName = ['ECG','Blood Pressure','Gastric Pressure','Breathing','EMG']
BChannel = 4;  
EKGChannel = 1;
GPChannel = 3;
BPChannel = 2;
laserChannel = 5;
%EMGChannel = 5


%% Analysis
ci = [1,2,7,8];
[fileDir] = uigetdir( 'Output');
cd(fileDir)
for iii = 1
    protocol_index=ci(iii);
    all_channel = find(protocol_set(ci(iii),FixInfo+1:end)>0)


    for jj = 1:(length(all_channel)-1)
        data_current = data_extract{protocol_index,all_channel(jj)};
        % Input parameter
        timeStart = 0
        timeEnd = protocol_set(protocol_index,4)*2.25;
        
        sampleRate = 1000/SampleRate;
        % sampleRate = sampleRate / binN;
        timePeriod = timeStart*sampleRate:timeEnd*sampleRate-1;
        T = 1/sampleRate;
        L = length(timePeriod);             % Length of signal
        t = (0:L-1)*T;
        nrep = 3

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

        windowOfInterest = 0.75
        lengthToAnalyze =  floor(laserDuration*windowOfInterest);
        analysisWindow = [];
        analysisWindow(:,1) =[ (laserDuration*2-lengthToAnalyze+1):1:(laserDuration*2)];  % periods during stimulation to be analyzed
        analysisWindow(:,2) =[ (laserDuration-lengthToAnalyze+1): 1:laserDuration];       % periods before stimulation to be analyzed

        BwindowOfInterest = 1
        BlengthToAnalyze =  floor(laserDuration*BwindowOfInterest);
        BanalysisWindow = [];
        BanalysisWindow(:,1) =[ (laserDuration+1):1:(laserDuration+BlengthToAnalyze)];  % periods during stimulation to be analyzed
        BanalysisWindow(:,2) =[ (laserDuration-BlengthToAnalyze+1): 1:laserDuration];       % periods before stimulation to be analyzed

        %
        
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

        %%%%%%%%%%  Breathing %%%%%%%%%%%%%
        if BChannel == all_channel(jj)
            % breathing parameters
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
            single_z_amp_tv ={};

            for k = 1:size(info,1)
                if k <0 %44
                    ifplot=1;
                else
                    ifplot=0;
                end
                for nn = 1:nrep
                    kk = (k-1)*nrep +nn;
                    maxBreathRate = Breathing_parameter{kk,1};  
                    %maxBreathRate = 6 % estimated maximal breath rate in Hz;
                    %pbreaththresholdlow = 0.20;      % recom: 0.18;
                    %pbreaththreshold = 0.8;         % recom: 0.5; (the lower the number, the higher the thd)
                    pbreaththresholdlow = Breathing_parameter{kk,3};      % recom: 0.18;
                    pbreaththreshold =Breathing_parameter{kk,2};
                    pSighThd = Breathing_parameter{kk,4};;% recom: 0.5; (the lower the number, the higher the thd)

                    breathData = detrend(data_current(:,kk));
                    for j = 1:2         % stim v. no stim
                        BreathingInWindow = breathData(BanalysisWindow(:,j));
                        [breathRate(kk,j),mVolume(kk,j),TVolume(kk,j),TVolume_auc(kk,j),TVolume_auc_nogasp(kk,j),mVolume_nogasp(kk,j),nSigh(kk,j), abspbthreshold,abspbthresholdlow,breathEvents] = analyzeBreath(ifplot,protocol_index, k, nn,BreathingInWindow,BanalysisWindow(:,j),sampleRate,maxBreathRate,pbreaththreshold,pbreaththresholdlow,pSighThd);     
                    end   
                    changeBreathRate(kk) = (breathRate(kk,1) - breathRate(kk,2))/(breathRate(kk,1) +breathRate(kk,2));
                    changeTVolume(kk) = (TVolume(kk,1) - TVolume(kk,2))/(TVolume(kk,1)+TVolume(kk,2));   
                        %     changeNSigh(i) = (nSigh(i,1) - nSigh(i,2))/nSigh(i,2) * 100;
                   
                    changeTVolume_auc(kk) = (TVolume_auc(kk,1) - TVolume_auc(kk,2))/(TVolume_auc(kk,1)+TVolume_auc(kk,2));   
                    changemVolume(kk) =(mVolume(kk,1) - mVolume(kk,2))/(mVolume(kk,1)+mVolume(kk,2));    
                    changeTVolume_auc_nogasp(kk) = (TVolume_auc_nogasp(kk,1) - TVolume_auc_nogasp(kk,2))/(TVolume_auc_nogasp(kk,1)+TVolume_auc_nogasp(kk,2));   
                    changemVolume_nogasp(kk) =(mVolume_nogasp(kk,1) - mVolume_nogasp(kk,2))/(mVolume_nogasp(kk,1)+mVolume_nogasp(kk,2));    
                    
                    changeNSigh(kk) = nSigh(kk,1) - nSigh(kk,2);

                    nnbin = 10; % Change!!! number of bins in each stage
                    bin = laserDuration/nnbin; % bin
                    tt = timeEnd*sampleRate;
                    bintime = bin/sampleRate;
                    %baseline = [1:nnbin];
                    baseline = [nnbin-9:nnbin];
                    
                    baseline_ms=BanalysisWindow(:,1);

                    [TV_auc_ci{kk},TV_auc_ci_nogasp{kk},single_z_amp_tv{kk},single_ci_amp_tv{kk},BR_CI{kk}, Delta_amp{kk}, CI_amp{kk}, Delta_amp_tv{kk}, CI_amp_tv{kk},BreathRate{kk},BR_zscore{kk},amp_stage{kk},amp_stage_tv{kk},amp_window{kk},amp_stage_raw{kk}] = analyzeBreath_SW(BlengthToAnalyze,protocol_index, k, nn,tt,bin,baseline,baseline_ms,laserDuration,breathData,timePeriod,sampleRate,maxBreathRate,pbreaththreshold,pbreaththresholdlow,pSighThd, abspbthreshold,abspbthresholdlow);
                    changeMinuteVolume(kk) = (sum(amp_window{kk}{:,1}) - sum(amp_window{kk}{:,2}))/(sum(amp_window{kk}{:,1})+sum(amp_window{kk}{:,2}));
                    single_z_amp_tv_trial{kk}= mean(single_z_amp_tv{kk});

                end
                for j =1:2
                    breathRate_mouse(k,j) = mean(breathRate((k-1)*3+1:(k-1)*3+nrep,j));
                    TVolume_mouse(k,j) = mean(TVolume((k-1)*3+1:(k-1)*3+nrep,j));
                    nSigh_mouse(k,j) = mean(nSigh((k-1)*3+1:(k-1)*3+nrep,j));

                    TVolume_auc_mouse(k,j) = mean(TVolume_auc((k-1)*3+1:(k-1)*3+nrep,j));
                    mVolume_mouse(k,j) = mean(mVolume((k-1)*3+1:(k-1)*3+nrep,j));
                    TVolume_auc_nogasp_mouse(k,j) = mean(TVolume_auc_nogasp((k-1)*3+1:(k-1)*3+nrep,j));
                    mVolume_nogasp_mouse(k,j) = mean(mVolume_nogasp((k-1)*3+1:(k-1)*3+nrep,j));
                end
                    single_z_amp_tv_mouse{k}= mean(cat(1,single_z_amp_tv{(k-1)*3+1:(k-1)*3+nrep}));
                
                    changeBreathRate_mouse(k) = mean(changeBreathRate((k-1)*3+1:(k-1)*3+nrep));
                    changeTVolume_mouse(k) = mean(changeTVolume((k-1)*3+1:(k-1)*3+nrep));
                    
                    changeTVolume_auc_mouse(k) = mean(changeTVolume_auc((k-1)*3+1:(k-1)*3+nrep));
                    changemVolume_mouse(k) = mean(changemVolume((k-1)*3+1:(k-1)*3+nrep));
                    changeTVolume_auc_nogasp_mouse(k) = mean(changeTVolume_auc_nogasp((k-1)*3+1:(k-1)*3+nrep));
                    changemVolume_nogasp_mouse(k) = mean(changemVolume_nogasp((k-1)*3+1:(k-1)*3+nrep));

                    changeNSigh_mouse(k) = mean(changeNSigh((k-1)* 3+1:(k-1)*3+nrep));
                    changeMinuteVolume_mouse(k) = mean(changeMinuteVolume((k-1)*3+1:(k-1)*3+nrep));
               
            end

           % write data
            clear writeBR; clear writeBR_mouse;
            writeBR = {'Trial #',...
                'Breath Rate Post Stim','Breath Rate Pre Stim','Breath Rate Change %',...
                'TV Post Stim','TV Pre Stim','TV Change %',...
                'TV AUC Post Stim','TV AUC Pre Stim','TV AUC Change %',...
                'MV Post Stim','MV Pre Stim','MV Change %',...
                'TV AUC nogasp Post Stim','TV AUC nogasp Pre Stim','TV AUC nogasp Change %',...
                'MV nogasp Post Stim','MV nogasp Pre Stim','MV nogasp Change %',...
                'nSigh Post Stim','nSigh Pre Stim','nSigh Change #',...
                'Minute Volume Change %'};
            for i = 1:size(info,1)*nrep
                writeBR = [writeBR;{i,...
                    breathRate(i,1),breathRate(i,2),changeBreathRate(i),...
                    TVolume(i,1),TVolume(i,2),changeTVolume(i),...
                    TVolume_auc(i,1),TVolume_auc(i,2),changeTVolume_auc(i),...
                    mVolume(i,1),mVolume(i,2),changemVolume(i),...
                    TVolume_auc_nogasp(i,1),TVolume_auc_nogasp(i,2),changeTVolume_auc_nogasp(i),...
                    mVolume_nogasp(i,1),mVolume_nogasp(i,2),changemVolume_nogasp(i),...
                    nSigh(i,1),nSigh(i,2),changeNSigh(i),changeMinuteVolume(i)}];
            end
            writeBR = [writeBR;{[], timeStart,timeEnd,sampleRate,...
                'pSighThd',pSighThd,[],...
                [],[],[],...
                [],[],[],...
                [],[],[],...
                [],[],[],...
                [],[],[],[]}]
            writeBR_mouse = {'Order','Mouse #','Group',...
                'Breath Rate Post Stim','Breath Rate Pre Stim','Breath Rate Change %',...
                'TV Post Stim','TV Pre Stim','TV Change %',...
                 'TV AUC Post Stim','TV AUC Pre Stim','TV AUC Change %',...
                'MV Post Stim','MV Pre Stim','MV Change %',...
                'TV AUC nogasp Post Stim','TV AUC nogasp Pre Stim','TV AUC nogasp Change %',...
                'MV nogasp Post Stim','MV nogasp Pre Stim','MV nogasp Change %',...
                'nSigh Post Stim','nSigh Pre Stim','nSigh Change #',...
                'Minute Volume Change %'};
           for i = 1:size(info,1)
                writeBR_mouse = [writeBR_mouse;{info{i,1},info{i,3},info{i,2},...
                    breathRate_mouse(i,1),breathRate_mouse(i,2),changeBreathRate_mouse(i),...
                    TVolume_mouse(i,1),TVolume_mouse(i,2),changeTVolume_mouse(i),...
                    TVolume_auc_mouse(i,1),TVolume_auc_mouse(i,2),changeTVolume_auc_mouse(i),...
                    mVolume_mouse(i,1),mVolume_mouse(i,2),changemVolume_mouse(i),...
                    TVolume_auc_nogasp_mouse(i,1),TVolume_auc_nogasp_mouse(i,2),changeTVolume_auc_nogasp_mouse(i),...
                    mVolume_nogasp_mouse(i,1),mVolume_nogasp_mouse(i,2),changemVolume_nogasp_mouse(i),...
                    nSigh_mouse(i,1),nSigh_mouse(i,2),changeNSigh_mouse(i),changeMinuteVolume_mouse(i)}];
            end
            writeBR_mouse = [writeBR_mouse;{[],[],[],...
                timeStart,timeEnd,sampleRate,...
                'pSighThd',pSighThd,[],...
                [],[],[],...
                [],[],[],...
                [],[],[],...
                [],[],[],...
                [],[],[],[]}];
            writecell(writeBR,['250809_Protocol',num2str(protocol_index),'_Breathing_trial.csv']);
            writecell(writeBR_mouse,['250809_Protocol',num2str(protocol_index),'_Breathing_mouse.csv']);

            matname =['250809_Protocol',num2str(protocol_index),'_Breathing_plot_data.mat']
            save (matname,'TV_auc_ci','TV_auc_ci_nogasp','single_z_amp_tv_mouse','single_z_amp_tv_trial','single_z_amp_tv','single_ci_amp_tv','BR_zscore','BR_CI','protocol_index', 'BreathRate', 'amp_stage','Delta_amp','CI_amp','amp_stage_tv','Delta_amp_tv','CI_amp_tv','amp_stage_raw','bin')
            bintime = bin/sampleRate;
        end
        % 
        % %%%%%%%%%  EKG %%%%%%%%%%%%%
        % if EKGChannel == all_channel(jj)
        %     maxHeartRate = 13;      % estimated maximal breath rate in Hz;
        %     EKGThd = -0.1;
        %     EKGThdLow = 0.1;
        % 
        %     HeartRate = {};
        %     amp_stage = {};
        %     HR_zscore = {};
        %     for k = 1: size(info,1)
        %         for nn = 1:nrep
        %             kk = (k-1)*nrep +nn;    
        % 
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
        %             heartRate_mouse(k,j) = mean(heartRate((k-1)*3+1:(k-1)*3+nrep,j));
        %             EKGAmplitude_mouse(k,j) = mean(EKGAmplitude((k-1)*3+1:(k-1)*3+nrep,j));
        %             changeHeartRate_mouse(k) = mean(changeHeartRate((k-1)*3+1:(k-1)*3+nrep));
        %             changeEKGAmplitude_mouse(k) = mean(changeEKGAmplitude((k-1)*3+1:(k-1)*3+nrep));
        %         end
        %     end
        % 
        %     % write data
        %     clear writeEKG; clear writeEKG_mouse;
        %     writeEKG = {'Trial #','Heart Rate Post Stim','Heart Rate Pre Stim','Heart Rate Change %','EKG Amp Post Stim','EKG Amp Pre Stim','EKG Amp Change %'};
        %     for i = 1:size(info,1)*nrep
        %         writeEKG = [writeEKG;{i,heartRate(i,1),heartRate(i,2),changeHeartRate(i),EKGAmplitude(i,1),EKGAmplitude(i,2),changeEKGAmplitude(i)}];
        %     end
        %     writeEKG = [writeEKG;{[],timeStart,timeEnd,sampleRate,'EKGThd',EKGThd,[]}];
        % 
        %     writeEKG_mouse = {'Mouse #','Heart Rate Post Stim','Heart Rate Pre Stim','Heart Rate Change %','EKG Amp Post Stim','EKG Amp Pre Stim','EKG Amp Change %'};
        %     for i = 1:size(info,1)
        %         writeEKG_mouse = [writeEKG_mouse;{i,heartRate_mouse(i,1),EKGAmplitude_mouse(i,2),changeHeartRate_mouse(i),EKGAmplitude_mouse(i,1),EKGAmplitude_mouse(i,2),changeEKGAmplitude_mouse(i)}];
        %     end
        %     writeEKG_mouse = [writeEKG_mouse;{[],timeStart,timeEnd,sampleRate,'EKGThd',EKGThd,[]}];
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

        % % 
        % %%%%%%%%%  Blood Pressure %%%%%%%%%%%%%
        % if BPChannel == all_channel(jj)
        %     EstMaxBP = 15 %Hz
        %     MinPeakH_max = 20
        %     MinPeakH_min = 150
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
        %     for k = 1: size(info,1)
        %         for nn = 1:nrep
        %             kk = (k-1)*nrep +nn
        %             % maxBreathRate = Breathing_parameter{kk,1};  
        %             % pbreaththresholdlow = Breathing_parameter{kk,2};      % recom: 0.18;
        %             % pbreaththreshold =Breathing_parameter{kk,3};
        % 
        %             BPData = data_current(:,kk);
        %             % figure
        %             % plot(BPData)
        %             % hold on
        %             % plot(bandstop(BPData,[0.2 3],sampleRate))
        %             BPData = bandstop(BPData,[0.2 3],sampleRate);
        %             for j = 1:2         % stim v. no stim
        %                 BPInWindow = BPData(BPanalysisWindow(:,j));
        %                 [BP_windowmean(kk,j),SBP_windowmean(kk,j),DBP_windowmean(kk,j),PulsePressure_windowmean(kk,j)] = analyzeBP(protocol_index, k, nn,BPInWindow,BPanalysisWindow(:,j),sampleRate,EstMaxBP,MinPeakH_max,MinPeakH_min);
        %             end   
        %             changeSBP(kk) = (SBP_windowmean(kk,1) - SBP_windowmean(kk,2))/(SBP_windowmean(kk,1) + SBP_windowmean(kk,2));
        %             changeDBP(kk) = (DBP_windowmean(kk,1) - DBP_windowmean(kk,2))/(DBP_windowmean(kk,1) + DBP_windowmean(kk,2));
        %             changePP(kk) = (PulsePressure_windowmean(kk,1) - PulsePressure_windowmean(kk,2))/(PulsePressure_windowmean(kk,1) + PulsePressure_windowmean(kk,2));
        %             deltaBP(kk)=(SBP_windowmean(kk,1)+DBP_windowmean(kk,1))/2-(SBP_windowmean(kk,2)+DBP_windowmean(kk,2))/2
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
        %             SBP_mouse(k,j) = mean(SBP_windowmean((k-1)*3+1:(k-1)*3+nrep,j));
        %             DBP_mouse(k,j) = mean(DBP_windowmean((k-1)*3+1:(k-1)*3+nrep,j));
        %             PulsePressure_mouse(k,j) = mean(PulsePressure_windowmean((k-1)*3+1:(k-1)*3+nrep,j));
        %         end
        %             changeSBP_mouse(k) = mean(changeSBP((k-1)*3+1:(k-1)*3+nrep));
        %             changeDBP_mouse(k) = mean(changeDBP((k-1)*3+1:(k-1)*3+nrep));
        %             changePP_mouse(k) = mean(changePP((k-1)*3+1:(k-1)*3+nrep));
        %             deltaBP_mouse(k) = mean(deltaBP((k-1)*3+1:(k-1)*3+nrep));
        %     end
        % 
        %    % write data
        %     clear writeBP; clear writeBP_mouse;
        %     writeBP = {'Trial #','deltaBP','SBP Post Stim','SBP Pre Stim','SBP Change %','DBP Post Stim','DBP Pre Stim','DBP Change %','PulsePresure Post Stim','PulsePresure Pre Stim','PulsePresure Change %'};
        %     for i = 1:size(info,1)*nrep
        %         writeBP = [writeBP;{i,deltaBP,SBP_windowmean(i,1),SBP_windowmean(i,2),changeSBP(i),DBP_windowmean(i,1),DBP_windowmean(i,2),changeDBP(i),PulsePressure_windowmean(i,1),PulsePressure_windowmean(i,2),changePP(i)}];
        %     end
        %     writeBP = [writeBP;{[],timeStart,timeEnd,sampleRate,[],[],[],[],[],[],[]}];
        % 
        %     writeBP_mouse = {'Mouse #','deltaBP','SBP Post Stim','SBP Pre Stim','SBP Change %','DBP Post Stim','DBP Pre Stim','DBP Change %','PulsePresure Post Stim','PulsePresure Pre Stim','PulsePresure Change %'};
        %     for i = 1:size(info,1)
        %         writeBP_mouse = [writeBP_mouse;{i, deltaBP_mouse(i),SBP_mouse(i,1),SBP_mouse(i,2),changeSBP_mouse(i),DBP_mouse(i,1),DBP_mouse(i,2),changeDBP_mouse(i),PulsePressure_mouse(i,1),PulsePressure_mouse(i,2),changePP_mouse(i)}];
        %     end
        %     writeBP_mouse = [writeBP_mouse;{[],timeStart,timeEnd,sampleRate,[],[],[],[],[],[],[]}];
        %     % [fileDir] = uigetdir( 'Output');
        %     % cd(fileDir);
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
        % 
        %     GP={}
        %     for k = 1: size(info,1)
        %         for nn = 1:nrep
        %             kk = (k-1)*nrep +nn
        %             % maxBreathRate = Breathing_parameter{kk,1};  
        %             % pbreaththresholdlow = Breathing_parameter{kk,2};      % recom: 0.18;
        %             % pbreaththreshold =Breathing_parameter{kk,3};
        % 
        %             samplepts = linspace(0,1,lengthToAnalyze);
        %             GPData = data_current(:,kk);
        %             GP{k,nn}= GPData
        %             if isempty(GPData)
        %                 comtinue
        %             end
        %             for j = 1:2         % stim v. no stim
        %                 GPInWindow = GPData(analysisWindow(:,j));  
        %                 [polorder,detrendGPdata(:,j),trendGP(:,j)] = RobustDetrend(GPInWindow,nDetrend,conflvl,samplepts);
        %                 detrendGPdata(:,j) = smooth(detrendGPdata(:,j),sampleRate);
        %                 %figure;plot(squeeze(detrendGPdata(:,j)));hold on;plot(squeeze(trendGP(:,j)));hold on; plot(GPInWindow)
        %                 %title(sprintf('Animal #%d Trial #%d, Stim #%d', k,nn,j));
        %                 phasic(kk,j) = sum(detrendGPdata(:,j)-min(detrendGPdata(:,j)));
        %                 tonic(kk,j) = mean(trendGP(:,j));
        %                 GPInWindow = smooth(GPInWindow);
        %                 % total(i,j) = sum(GPInWindow-min(GPInWindow));
        %                 total(kk,j) = sum(GPInWindow);
        %                 total_mean(kk,j)=mean(GPInWindow);
        %             end  
        %             changePhasic(kk) = (phasic(kk,1) - phasic(kk,2))/(phasic(kk,1) + phasic(kk,2));
        %             changeTonic(kk) = (tonic(kk,1) - tonic(kk,2))/(tonic(kk,1) + tonic(kk,2));
        %             changeTotal(kk) = (total(kk,1) - total(kk,2))/(total(kk,1) + total(kk,2));
        %             changeTotal_mean(kk) = (total_mean(kk,1) - total_mean(kk,2))/ (total_mean(kk,1) + total_mean(kk,2));
        %         end
        %         for j =1:2
        %             phasic_mouse(k,j) = mean(phasic((k-1)*3+1:(k-1)*3+nrep,j));
        %             tonic_mouse(k,j) = mean(tonic((k-1)*3+1:(k-1)*3+nrep,j));
        %             total_mouse(k,j) = mean(total((k-1)*3+1:(k-1)*3+nrep,j));
        %             total_mean_mouse(k,j) = mean(total_mean((k-1)*3+1:(k-1)*3+nrep,j));
        %         end
        % 
        %         changePhasic_mouse(k) = mean(changePhasic((k-1)*3+1:(k-1)*3+nrep));
        %         changeTonic_mouse(k) = mean(changeTonic((k-1)*3+1:(k-1)*3+nrep));
        %         changeTotal_mouse(k) = mean(changeTotal((k-1)*3+1:(k-1)*3+nrep));
        %         changeTotal_mean_mouse(k) = mean(changeTotal_mean((k-1)*3+1:(k-1)*3+nrep));
        % 
        %     end
        %     % write data
        %     clear writeGP;clear writeGP_mouse;
        %     writeGP = {'Trial #','Phsic Post Stim','Phasic Pre Stim','Phsic Change %','Tonic Post Stim','Tonic Pre Stim','Tonic Change %','Total Post Stim','Total Pre Stim','Total Change %'};
        %     for i = 1:size(info,1)*nrep
        %         writeGP = [writeGP;{i,phasic(i,1),phasic(i,2),changePhasic(i),tonic(i,1),tonic(i,2),changeTonic(i),total(i,1),total(i,2),changeTotal(i)}];
        %     end
        %     writeGP = [writeGP;{[],timeStart,timeEnd,sampleRate,[],[],[],[],[],[]}];
        % 
        %     writeGP_mouse = {'Mouse #','Phsic Post Stim','Phasic Pre Stim','Phsic Change %','Tonic Post Stim','Tonic Pre Stim','Tonic Change %','Total Post Stim','Total Pre Stim','Total Change %'};
        %     for i = 1:size(info,1)
        %         writeGP_mouse = [writeGP_mouse;{i,phasic_mouse(i,1),phasic_mouse(i,2),changePhasic_mouse(i),tonic_mouse(i,1),tonic_mouse(i,2),changeTonic_mouse(i),total_mouse(i,1),total_mouse(i,2),changeTotal_mouse(i)}];
        %     end
        %     writeGP_mouse = [writeGP_mouse;{[],timeStart,timeEnd,sampleRate,[],[],[],[],[],[]}];
        %     % [fileDir] = uigetdir( 'Output');
        %     % cd(fileDir);
        %     writecell(writeGP,['Protocol',num2str(protocol_index),'_GastricPressure_trial.csv']);
        %     writecell(writeGP_mouse,['Protocol',num2str(protocol_index),'_GastricPressure_mouse.csv']);
        % 
        %     matname =['Protocol',num2str(protocol_index),'_GastricPressure_plot_data.mat'];
        %     save (matname,'GP','protocol_index');
        %     %bintime = bin/sampleRate;
        % end

    end
end




            


   