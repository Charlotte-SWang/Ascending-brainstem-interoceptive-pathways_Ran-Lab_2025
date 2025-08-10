function [heartRate,EKGAmplitude] = analyzeEKG(protocol_index, k,EKGData,timePeriod,sampleRate,maxHeartRate,EKGThd,EKGThdLow)
% Ran, analysis of breath rate and tidal volume
    minPeakDis = 1/maxHeartRate*sampleRate;
    [VEKG,EKGevents] = findpeaks(EKGData,timePeriod,'MinPeakDistance',minPeakDis,'MinPeakHeight',EKGThd);
    invpresData = -EKGData;
    [junk, lk5] = findpeaks(invpresData,timePeriod,'MinPeakDistance',minPeakDis,'MinPeakHeight',-EKGThdLow);
    temp = [];

    for iInd = 1:length(lk5)   
        temp(iInd) = find(timePeriod == lk5(iInd),1,'first');
    end
    if isempty(temp) == 1
        VEKGLow = 0;
    else
        VEKGLow = EKGData(temp);
        %if k==29 | k ==25|k== 24|k== 15
            figure;
            plot(timePeriod,EKGData,EKGevents,VEKG,'ro');
            hold on;
            plot(timePeriod,EKGData,lk5,VEKGLow,'r*');
            title(['protocol', num2str(protocol_index), ' mouse',num2str(k),'EKG analysis']);
        %end;

        % use the following lines when plotting a part of the plot
        timeLength = length(timePeriod);
        timeLimit = floor(timeLength/100);
        tempEKGevents = EKGevents(EKGevents<timePeriod(timeLimit));
        templk5 = lk5(lk5<timePeriod(timeLimit));
        tempVEKG = VEKG(1:length(tempEKGevents));
        tempVEKGLow = VEKGLow(1:length(templk5));
        % figure;
        % plot(timePeriod(1:timeLimit),EKGData(1:timeLimit),tempEKGevents,tempVEKG,'ro');
        % hold on;
        % plot(timePeriod(1:timeLimit),EKGData(1:timeLimit),templk5,tempVEKGLow,'r*');
        % title('EKG analysis');
    % end
    aveVEKG = mean(VEKG);
    aveVEKGLow = mean(VEKGLow);
    EKGAmplitude = aveVEKG - aveVEKGLow;
    heartRate = length(VEKGLow)/(length(timePeriod) / sampleRate);
end

