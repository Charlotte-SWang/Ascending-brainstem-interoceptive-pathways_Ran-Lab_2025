function [HeartRate,HR_zscore,amp_stage,HR_CI] = analyzeEKG_SW(protocol_index, k, nn,tt,bin,baseline,laserDuration,EKGData,timePeriod,sampleRate,maxHeartRate,EKGThd,EKGThdLow)
% Ran, analysis of breath rate and tidal volume
    minPeakDis = 1/maxHeartRate*sampleRate;
    [VEKG,EKGevents] = findpeaks(EKGData,timePeriod,'MinPeakDistance',minPeakDis,'MinPeakHeight',EKGThd);
    invpresData = -EKGData;
    locs_max = EKGevents;
    [junk, locs_min] = findpeaks(invpresData,timePeriod,'MinPeakDistance',minPeakDis,'MinPeakHeight',-EKGThdLow);
    pks_max = VEKG;
    pks_min =junk;

    nbin = ceil(tt/bin);
    numberformin = 60*sampleRate/bin; %s to ms 
    HR = histcounts(locs_max,'BinWidth',bin)*numberformin;
    HeartRate = HR;
    z_m = mean(HR(baseline));
    z_sd = std(HR(baseline));
    HR_zscore = (HR-z_m)./z_sd;
    HR_CI = (HR-z_m)./(HR+z_m);
    la=length(locs_max);
    lb=length(locs_min);
    % Amplitude analysis
    while la ~= lb
        if length(locs_max) > length(locs_min)
           locs_min(length(locs_min)+1:length(locs_max))= [inf]*(length(locs_max) - length(locs_min));
           locs_diff =locs_min-locs_max;
           locs_max(find(locs_diff>0,1,'first'))= [];
           pks_max(find(locs_diff>0,1,'first'))= [];
           locs_min(find(locs_min==inf))=[];
           la=length(locs_max);
           lb=length(locs_min);
        else
           locs_max(length(locs_max)+1:length(locs_min))= [-inf]*(length(locs_min) - length(locs_max));
           locs_diff =locs_min-locs_max;
           locs_min(find(locs_diff>0,1,'first'))= [];
           pks_min(find(locs_diff>0,1,'first'))= [];
           locs_max(find(locs_max==-inf))=[];
           la=length(locs_max);
           lb=length(locs_min);
        end
    end

    peak_maxmin = [pks_max pks_min];
    amp = abs(sum(peak_maxmin,2));
    amp_idx = locs_min;
    amp_stage = [];

    amp_stage{1,1} = amp(find(amp_idx <= laserDuration));
        %On
    amp_stage{1,2} = amp(find(amp_idx <= laserDuration*2 & amp_idx > (laserDuration)));
    amp_stage{1,3} = amp(find(amp_idx >laserDuration*2));

    
end

