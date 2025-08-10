%% analyze physiology of the PVH and VLM stim experiment
% Ran 20201019
% Adapted by Shiqi 20250116
clear;close all;clc;
addpath('E:\tempAnalysis\matlabCodesScripps\functions');
%% Input parameters
timeStart = 18*60;
% timeEnd = timeStart + 500;
timeEnd = 36*60;
sampleRate = 1000;
fileDir = ('E:\tempAnalysis\2023opto\');
animalID = '20 VLME';
 
tempAnalysisCutoff = 2000; % (when starting time is >2000, it is blood pressure that is being analyzed).
fiveSecBreathingAnalysiCutoff = 200; % (when starting time is <200, it is 5-sec breathing that is being analyzed).

% channels
% Modified for Small Volume Injection 2025/01
ChannelName = ['ECG','Blood Pressure','Gastric Pressure','Breathing','EMG']
BChannel = 4;  
EKGChannel = 1;
GPChannel = 3;
BPChannel = 2;
laserChannel = 6;
EMGChannel = 5

% EMGChannel = 4;
% BChannel = 2;
% % EKGChannel = 2;
% % GPChannel = 1;
% % BPChannel = 3;
% laserChannel = 4;
% % EMGChannel = 4;

% breathing parameters
maxBreathRate = 6;      % estimated maximal breath rate in Hz;
% pbreaththresholdlow = 0.32;      % recom: 0.18;
% pbreaththreshold = 0.38;         % recom: 0.5; (the lower the number, the higher the thd)
pbreaththresholdlow = 0.18;      % recom: 0.18;
pbreaththreshold = 0.5;         % recom: 0.5; (the lower the number, the higher the thd)
pSighThd = 4;                  % recom: 0.5

% EKG parameters
maxHeartRate = 13;      % estimated maximal breath rate in Hz;
EKGThd = 0.1;
EKGThdLow = -0.3;

% laser parameters
if timeStart < fiveSecBreathingAnalysiCutoff        % when analyzing 5 seconds breathing or 
    minLaserDistance = 6;                           % a number higher than the estimated laser stimulation train band width in second;
elseif timeStart > tempAnalysisCutoff               % use these five lines when analyzing blood pressure. Use the 6th line when analyzing others.
    minLaserDistance = 70;
else
    % minLaserDistance = 6;
    minLaserDistance = 310;
end


% GP parameters
lowpassCutoffGP = 1.5;
arbitraryBcgdPeriod = 5;
nDetrend = 3;
conflvl = 0.9;



fileDir = [fileDir animalID];
fileName = dir([fileDir '\*.mat']);
fileDirName = [fileDir '\' fileName(1).name];
load(fileDirName);
% sampleRate = sampleRate / binN;
timePeriod = timeStart*sampleRate:timeEnd*sampleRate;
T = 1/sampleRate;
L = length(timePeriod);             % Length of signal
t = (0:L-1)*T;

[pathname1] = uigetdir( 'Select the physiology data');
cd(pathname1);

%% laser onset and offset

laserData = data(timePeriod,laserChannel);
minLaserDistance = minLaserDistance * sampleRate;
laserData(laserData >= 4) = 4;
[plaseron,onsets] = findpeaks(laserData,timePeriod,'MinPeakDistance',minLaserDistance,'MinPeakHeight',2);
if length(onsets) ~= 1                  % only one trial
    laserDuration = (onsets(2) - onsets(1))/2;
else
    error('manual input of laser duration is required');
% laserDuration = 300244;                     % Laser duration for 5 min stim
end
laserDuration = floor(laserDuration);

offsets = onsets + laserDuration;
figure;
plot(timePeriod,laserData,onsets,plaseron,'ro',offsets,plaseron,'ko');
% close all;

if length(onsets) > 10
    error('Check minLaserDistance, sampleRate, and timePeriod!');
end

minLength = min(offsets-onsets);        % make sure the time periods analyzed has equal length
lengthToAnalyze = floor(minLength*windowOfInterest); % length of the time periods to be analyzed
analysisWindow = zeros(length(onsets),lengthToAnalyze,2);
nTrial = length(min(onsets,offsets)); % number of trials to analyze
for i = 1:nTrial
    analysisWindow(i,:,1) = [offsets(i) - lengthToAnalyze+1:offsets(i)];    % periods during stimulation to be analyzed
    analysisWindow(i,:,2) = [onsets(i) - lengthToAnalyze+1:onsets(i)];      % periods before stimulation to be analyzed
end
analysisWindow = analysisWindow - timePeriod(1);






if timeStart < fiveSecBreathingAnalysiCutoff
    % breathing
    breathData = -data(timePeriod,BChannel);
    for j = 1:2         % stim v. no stim
        for i = 1:nTrial
            BreathingInWindow = breathData(analysisWindow(i,:,j));
            [breathRate(i,j),TVolume(i,j), nSigh(i,j), abspbthreshold,abspbthresholdlow,breathEvents] = analyzeBreath(BreathingInWindow,analysisWindow(i,:,j),sampleRate,maxBreathRate,pbreaththreshold,pbreaththresholdlow,pSighThd);
        end
        breathRate(i+1,j) = mean(breathRate(1:nTrial,j));
        TVolume(i+1,j) = mean(TVolume(1:nTrial,j));
        nSigh(i+1,j) = mean(nSigh(1:nTrial,j));
    end
    for i = 1:nTrial
        changeBreathRate(i) = (breathRate(i,1) - breathRate(i,2))/breathRate(i,2) * 100;
        changeTVolume(i) = (TVolume(i,1) - TVolume(i,2))/TVolume(i,2) * 100;
        changeNSigh(i) = (nSigh(i,1) - nSigh(i,2))/nSigh(i,2) * 100;
        changeNSigh(i) = nSigh(i,1) - nSigh(i,2);
    end
    changeBreathRate(i+1) = mean(changeBreathRate);
    changeTVolume(i+1) = mean(changeTVolume);
    changeNSigh(i+1) = mean(changeNSigh);
    % write data
    clear writeBR;
    writeBR = {'Trial #','Breath Rate Post Stim','Breath Rate Pre Stim','Breath Rate Change %','TV Post Stim','TV Pre Stim','TV Change %','nSigh Post Stim','nSigh Pre Stim','nSigh Change #'};
    for i = 1:nTrial
        writeBR = [writeBR;{i,breathRate(i,1),breathRate(i,2),changeBreathRate(i),TVolume(i,1),TVolume(i,2),changeTVolume(i),nSigh(i,1),nSigh(i,2),changeNSigh(i)}];
    end
    writeBR = [writeBR;{'Mean',breathRate(i+1,1),breathRate(i+1,2),changeBreathRate(i+1),TVolume(i+1,1),TVolume(i+1,2),changeTVolume(i+1),nSigh(i+1,1),nSigh(i+1,2),changeNSigh(i+1)}];
    writeBR = [writeBR;{animalID,timeStart,timeEnd,sampleRate,'pSighThd',pSighThd,[],[],[],[]}];
    cd(fileDir);
    writecell(writeBR,['BR' animalID '.csv']);






elseif timeStart > tempAnalysisCutoff
    %% Blood pressure
    BPData = data(timePeriod,BPChannel);
    for j = 1:2         % stim v. no stim
        for i = 1:nTrial
            maxBP(i,j) = max(smooth(BPData(analysisWindow(i,:,j))));
        end
        maxBP(i+1,j) = mean(maxBP(1:nTrial,j));
    end
    for i = 1:nTrial
        changeBP(i) = (maxBP(i,1) - maxBP(i,2))/maxBP(i,2) * 100;
    end
    changeBP(i+1) = mean(changeBP);
    % write data
    writeBP = {'Trial #','Blood Pressure Post Stim','Blood Pressure Pre Stim','Blood Pressure Change %'};
    for i = 1:nTrial
        writeBP = [writeBP;{i,maxBP(i,1),maxBP(i,2),changeBP(i)}];
    end
    writeBP = [writeBP;{i,maxBP(i+1,1),maxBP(i+1,2),changeBP(i+1)}];
    writeBP = [writeBP;{animalID,timeStart,timeEnd,sampleRate}];
    cd(fileDir);
    writecell(writeBP,['BP' animalID '.csv']);


else
    %% Gastric pressure
    GPData = data(timePeriod,GPChannel);
    samplepts = linspace(0,1,lengthToAnalyze);
    for j = 1:2         % stim v. no stim
        for i = 1:nTrial
            GPInWindow = GPData(analysisWindow(i,:,j));
            [polorder,detrendGPdata(i,:,j),trendGP(i,:,j)] = RobustDetrend(GPInWindow,nDetrend,conflvl,samplepts);
            detrendGPdata(i,:,j) = smooth(detrendGPdata(i,:,j),sampleRate);
            figure;plot(squeeze(detrendGPdata(i,:,j)));hold on;plot(squeeze(trendGP(i,:,j)));
            title(sprintf('Trial #%d, Stim #%d', i,j));
            phasic(i,j) = sum(detrendGPdata(i,:,j)-min(detrendGPdata(i,:,j)));
            tonic(i,j) = mean(trendGP(i,:,j));
            GPInWindow = smooth(GPInWindow);
            % total(i,j) = sum(GPInWindow-min(GPInWindow));
            total(i,j) = sum(GPInWindow);
        end
        phasic(i+1,j) = mean(phasic(1:nTrial,j));
        tonic(i+1,j) = mean(tonic(1:nTrial,j));
        total(i+1,j) = mean(total(1:nTrial,j));
    end
    for i = 1:nTrial
        changePhasic(i) = (phasic(i,1) - phasic(i,2))/phasic(i,2) * 100;
        changeTonic(i) = (tonic(i,1) - tonic(i,2))/tonic(i,2) * 100;
        changeTotal(i) = (total(i,1) - total(i,2))/total(i,2) * 100;
    end
    changePhasic(i+1) = mean(changePhasic);
    changeTonic(i+1) = mean(changeTonic);
    changeTotal(i+1) = mean(changeTotal);
    % write data
    clear writeGP;
    writeGP = {'Trial #','Phsic Post Stim','Phasic Pre Stim','Phsic Change %','Tonic Post Stim','Tonic Pre Stim','Tonic Change %','Total Post Stim','Total Pre Stim','Total Change %'};
    for i = 1:nTrial
        writeGP = [writeGP;{i,phasic(i,1),phasic(i,2),changePhasic(i),tonic(i,1),tonic(i,2),changeTonic(i),total(i,1),total(i,2),changeTotal(i)}];
    end
    writeGP = [writeGP;{'Mean',phasic(i+1,1),phasic(i+1,2),changePhasic(i+1),tonic(i+1,1),tonic(i+1,2),changeTonic(i+1),total(i+1,1),total(i+1,2),changeTotal(i+1)}];
    writeGP = [writeGP;{animalID,timeStart,timeEnd,sampleRate,[],[],[],[],[],[]}];
    cd(fileDir);
    writecell(writeGP,[animalID 'GP.csv']);


    %% EKG
    EKGData = data(timePeriod,EKGChannel);
    minPeakDis = 1/maxBreathRate*sampleRate;
    for j = 1:2         % stim v. no stim
        for i = 1:nTrial
            EKGInWindow = EKGData(analysisWindow(i,:,j));
            [heartRate(i,j),EKGAmplitude(i,j)] = analyzeEKG(EKGInWindow,analysisWindow(i,:,j),sampleRate,maxHeartRate,EKGThd,EKGThdLow);
        end
        heartRate(i+1,j) = mean(heartRate(1:nTrial,j));
        EKGAmplitude(i+1,j) = mean(EKGAmplitude(1:nTrial,j));
    end
    for i = 1:nTrial
        changeHeartRate(i) = (heartRate(i,1) - heartRate(i,2))/heartRate(i,2) * 100;
        changeEKGAmplitude(i) = (EKGAmplitude(i,1) - EKGAmplitude(i,2))/EKGAmplitude(i,2) * 100;
    end
    changeHeartRate(i+1) = mean(changeHeartRate);
    changeEKGAmplitude(i+1) = mean(changeEKGAmplitude);

    % write data
    writeEKG = {'Trial #','Heart Rate Post Stim','Heart Rate Pre Stim','Heart Rate Change %','EKG Amp Post Stim','EKG Amp Pre Stim','EKG Amp Change %'};
    for i = 1:nTrial
        writeEKG = [writeEKG;{i,heartRate(i,1),heartRate(i,2),changeHeartRate(i),EKGAmplitude(i,1),EKGAmplitude(i,2),changeEKGAmplitude(i)}];
    end
    writeEKG = [writeEKG;{'Mean',heartRate(i+1,1),heartRate(i+1,2),changeHeartRate(i+1),EKGAmplitude(i+1,1),EKGAmplitude(i+1,2),changeEKGAmplitude(i+1)}];
    writeEKG = [writeEKG;{animalID,timeStart,timeEnd,sampleRate,'pSighThd',pSighThd,[]}];
    cd(fileDir);
    writecell(writeEKG,['EK' animalID '.csv']);

    %% breathing
    breathData = -data(timePeriod,BChannel);
    for j = 1:2         % stim v. no stim
        for i = 1:nTrial
            BreathingInWindow = breathData(analysisWindow(i,:,j));
            [breathRate(i,j),TVolume(i,j), nSigh(i,j), abspbthreshold,abspbthresholdlow,breathEvents] = analyzeBreath(BreathingInWindow,analysisWindow(i,:,j),sampleRate,maxBreathRate,pbreaththreshold,pbreaththresholdlow,pSighThd);
        end
        breathRate(i+1,j) = mean(breathRate(1:nTrial,j));
        TVolume(i+1,j) = mean(TVolume(1:nTrial,j));
        nSigh(i+1,j) = mean(nSigh(1:nTrial,j));
    end
    for i = 1:nTrial
        changeBreathRate(i) = (breathRate(i,1) - breathRate(i,2))/breathRate(i,2) * 100;
        changeTVolume(i) = (TVolume(i,1) - TVolume(i,2))/TVolume(i,2) * 100;
    %     changeNSigh(i) = (nSigh(i,1) - nSigh(i,2))/nSigh(i,2) * 100;
        changeNSigh(i) = nSigh(i,1) - nSigh(i,2);
    end
    changeBreathRate(i+1) = mean(changeBreathRate);
    changeTVolume(i+1) = mean(changeTVolume);
    changeNSigh(i+1) = mean(changeNSigh);
    % write data
    clear writeBR;
    writeBR = {'Trial #','Breath Rate Post Stim','Breath Rate Pre Stim','Breath Rate Change %','TV Post Stim','TV Pre Stim','TV Change %','nSigh Post Stim','nSigh Pre Stim','nSigh Change #'};
    for i = 1:nTrial
        writeBR = [writeBR;{i,breathRate(i,1),breathRate(i,2),changeBreathRate(i),TVolume(i,1),TVolume(i,2),changeTVolume(i),nSigh(i,1),nSigh(i,2),changeNSigh(i)}];
    end
    writeBR = [writeBR;{'Mean',breathRate(i+1,1),breathRate(i+1,2),changeBreathRate(i+1),TVolume(i+1,1),TVolume(i+1,2),changeTVolume(i+1),nSigh(i+1,1),nSigh(i+1,2),changeNSigh(i+1)}];
    writeBR = [writeBR;{animalID,timeStart,timeEnd,sampleRate,'pSighThd',pSighThd,[],[],[],[]}];
    cd(fileDir);
    writecell(writeBR,['B5' animalID '.csv']);
end
