function [ TV_auc_ci,TV_auc_ci_nogasp,single_z_amp_tv,single_ci_amp_tv,BR_CI, Delta_amp, CI_amp, Delta_amp_tv, CI_amp_tv, BreathRate,BR_zscore,amp_stage,amp_stage_tv,amp_window,amp_stage_raw] = analyzeBreath_SW(BlengthToAnalyze,protocol_index, k, nn,tt,bin,baseline,baseline_ms,laserDuration,breathData,timePeriod,sampleRate,maxBreathRate,pbreaththreshold,pbreaththresholdlow,pSighThd, abspbthreshold,abspbthresholdlow)
% Ran, analysis of breath rate and tidal volume
% Modified by Shiqi Wang
minPeakDis = 1/maxBreathRate*sampleRate;
abspbreaththreshold = pbreaththreshold;
abspbreaththresholdlow = pbreaththresholdlow;

% defining minPeakHeight
[junkpbreath,lk] = findpeaks(breathData,timePeriod,'MinPeakDistance',minPeakDis,'MinPeakHeight',abspbreaththreshold);

invpresData = -breathData;
[junk, lk2] = findpeaks(invpresData,timePeriod,'MinPeakDistance',minPeakDis,'MinPeakHeight',-abspbreaththresholdlow);
lk3 = [];
for iInd = 1:length(lk2)   
    lk3(iInd) = find(timePeriod == lk2(iInd));
end;
junkpbreathlow = breathData(lk3);
avepbreath = mean(junkpbreath);
avepbreathlow = mean(junkpbreathlow);
tidalVolume = avepbreath - avepbreathlow;
    
%% analyzing breath rate and tidal volume

% sigh
sighThd = pSighThd;
%sighThd = avepbreathlow + tidalVolume * pSighThd;
[psigh,lk4] = findpeaks(breathData,timePeriod,'MinPeakDistance',minPeakDis,'MinPeakHeight',sighThd);
% figure;
% plot(timePeriod,breathData,lk4,psigh,'ko');
% line([timePeriod(1) timePeriod(end)],[sighThd sighThd],'color',[0 0 0],'linestyle',':');
% %title('sigh analysis');
% title(['protocol',num2str(protocol_index),'mouse',num2str(k),' trial',num2str(nn),' Sigh analysis']);

nSigh = length(lk4);
abspbthd = abspbreaththreshold;
abspbthdlow = abspbreaththresholdlow;

% breath

[pks_max,locs_max] = findpeaks(breathData,timePeriod,'MinPeakDistance',minPeakDis,'MinPeakHeight',abspbreaththreshold);
    
invpresData = -breathData;
[pks_min,locs_min] = findpeaks(invpresData,timePeriod,'MinPeakDistance',minPeakDis,'MinPeakHeight',-abspbreaththresholdlow);

        %BreathStamp_max = Breath(BaselineStart:OffEnd).*0;
        %BreathStamp_max(locs_max,1) = 1;
         % ms
nbin = ceil(tt/bin);
numberformin = 60*sampleRate/bin; %s to ms 
BR = histcounts(locs_max,nbin)*numberformin;
BreathRate = BR;

z_m = mean(BR(baseline));
z_sd = std(BR(baseline));

BR_zscore = (BR-z_m)./z_sd;
BR_CI = (BR-z_m)./(BR+z_m);



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
    gasp=[];
    nogasp_loc = [];
    gasp_loc = [];
    gasp_total=[];
    sigh_idx = [];
    gasp_index_high=locs_max(find(pks_max>pSighThd)); 

    closest_smaller = nan(size( gasp_index_high)); % preallocate with NaN (if none found)

    for k = 1:length (gasp_index_high)
        smaller_vals = locs_min(locs_min <  gasp_index_high(k));  % only smaller values
        if ~isempty(smaller_vals)
            closest_smaller(k) = max(smaller_vals); % pick the largest of the smaller ones
        end
    end
    gasp_index_low=unique(closest_smaller);

    sigh_idx = sort([gasp_index_high gasp_index_low]);
    pks_min(pks_max>pSighThd)=[];
    locs_min(pks_max>pSighThd)=[];
    locs_max(pks_max>pSighThd)=[];
    pks_max(pks_max>pSighThd)=[];

    peak_maxmin = [pks_max pks_min];

    amp = abs(pks_max+pks_min);
    amp_idx = locs_max;

    % gasp(1)=0;
    % gasp(length(amp))=0;
    % for ii =2:length(amp)-1
    %     if amp(ii)>(amp(ii+1)*1.5) & amp(ii)>(amp(ii-1)*1.5)
    %         gasp(ii)=1;
    %     else
    %         gasp(ii)=0;
    %     end
    % end
    % 
    % amp_filter= abs(sum(peak_maxmin_filter,2));
    %     amp_idx_filter = locs_min(nogasp_loc);
 
    %%% amp_volume
    amp_stage = [];
    Delta_amp = [];
    CI_amp = [];
    for jj = 1:nbin
        return_idx = find(amp_idx <= jj*bin & amp_idx > (jj-1)*bin);
        if isempty(return_idx)
            amp_stage_raw(jj) =0;
        else
            amp_stage_raw(jj) = sum(amp(return_idx));
        end
    end
    z_m = mean(amp_stage_raw(baseline));
    z_sd = std(amp_stage_raw(baseline));
    amp_stage = (amp_stage_raw-z_m)./z_sd;
    Delta_amp = amp_stage_raw-z_m ;
    CI_amp = (amp_stage_raw-z_m )./(amp_stage_raw+z_m );
    %%% tidal_volume

    amp_stage_tv = [];
    Delta_amp_tv = [];
    CI_amp_tv = [];
    for jj = 1:nbin
        return_idx = find(amp_idx <= jj*bin & amp_idx > (jj-1)*bin);
        if isempty(return_idx)
            amp_stage_raw(jj) =NaN;
        else
            amp_stage_raw(jj) = mean(nonzeros(amp(return_idx)));
        end
        %amp_stage_raw(jj) = mean(asmp(find(amp_idx <= jj*bin & amp_idx > (jj-1)*bin)));
    end
    z_m = mean(amp_stage_raw(baseline),"omitnan");
    z_sd = std(amp_stage_raw(baseline),"omitnan");
    amp_stage_tv = (amp_stage_raw-z_m)./z_sd;
    Delta_amp_tv = amp_stage_raw-z_m ;
    CI_amp_tv = (amp_stage_raw-z_m )./(amp_stage_raw+z_m );

    single_amp_baseline=0;
    return_idx_tv = find(amp_idx >= (baseline(1)-1)*bin & amp_idx <= (baseline(end)*bin));
    single_amp_baseline = mean(nonzeros(amp(return_idx_tv)));
    return_stimulation_tv = find(amp_idx > baseline(end)*bin & amp_idx <= ((2*baseline(end)-baseline(1)+1)*bin)); 
    %amp_stage_raw(jj) = mean(asmp(find(amp_idx <= jj*bin & amp_idx > (jj-1)*bin)));
    single_ci_amp_tv=[];
    single_z_amp_tv=[];
    single_ci_amp_tv=(amp(return_stimulation_tv)-single_amp_baseline)./(amp(return_stimulation_tv)+single_amp_baseline);
    zz_sd = std(nonzeros(amp(return_idx_tv)),"omitnan");
    single_z_amp_tv=(amp(return_stimulation_tv)-single_amp_baseline)./zz_sd;

    amp_window = [];
    amp_window{1,1} = amp(find(amp_idx > laserDuration+1 & amp_idx <= laserDuration+BlengthToAnalyze));
    amp_window{1,2} = amp(find(amp_idx > laserDuration-BlengthToAnalyze & amp_idx <= laserDuration));
    
    %%% Calculate data for area under curve
    bD=breathData(baseline_ms);
    TV_auc=[];
    auc=[];
    bD_max = max(bD);
    bD_min = min(bD);
    zc_indices = find(bD(1:end-1) .* bD(2:end) < 0);
    for i = 1:length(zc_indices)-1
        idx_start = zc_indices(i);
        idx_end = zc_indices(i+1);
        auc(i) = sum(bD(idx_start:idx_end));
        % Breath volume change
    end
    
    if bD_max>abs(bD_min)
        if auc(1)<0
            j=1;
        else
            j=2;
        end
    else
        if auc(1)>0
            j=1;
        else
            j=2;
        end
    end
    while j+1<=length(auc)
        TV_auc(end+1)=abs(auc(j))+abs(auc(j+1));
        j=j+2;
    end
    TV_auc_baseline=mean(TV_auc);


    bD=breathData(baseline_ms(end)+1:end);
    TV_auc=[];
    auc=[];
    bD_max = max(bD);
    bD_min = min(bD);
    zc_indices = find(bD(1:end-1) .* bD(2:end) < 0);
    for i = 1:length(zc_indices)-1
        idx_start = zc_indices(i);
        idx_end = zc_indices(i+1);
        auc(i) = sum(bD(idx_start:idx_end));
        % Breath volume change
    end
    
    if bD_max>abs(bD_min)
        if auc(1)<0
            j=1;
        else
            j=2;
        end
    else
        if auc(1)>0
            j=1;
        else
            j=2;
        end
    end
    while j+1<=length(auc)
        TV_auc(end+1)=abs(auc(j))+abs(auc(j+1));
        j=j+2;
    end
    TV_auc_ci=(TV_auc-TV_auc_baseline)./(TV_auc+TV_auc_baseline);
%%%%%% nogasp
    checksigh=1
    %%% Exclude gasp
    bD=breathData(baseline_ms);
    TV_auc_nogasp=[];
    auc_nogasp=[];
    bD_max = max(bD);
    bD_min = min(bD);
    if length(sigh_idx)>0
        for i = 1:length(zc_indices)-1
            idx_start = zc_indices(i);
            idx_end = zc_indices(i+1);
            if idx_end<sigh_idx(checksigh)
                auc_nogasp(i) = sum(bD(idx_start:idx_end));
            elseif idx_start>sigh_idx(checksigh)
                auc_nogasp(i) = sum(bD(idx_start:idx_end));
            elseif (idx_end>sigh_idx(checksigh)) & (idx_start<sigh_idx(checksigh)) %exclude sigh
                if checksigh < length(sigh_idx)
                    checksigh = checksigh + 1;
                end
                continue
            end
            % Breath volume change
        end
    else
       
        for i = 1:length(zc_indices)-1
            idx_start = zc_indices(i);
             idx_end = zc_indices(i+1);
             auc_nogasp(i) = sum(breathData(idx_start:idx_end));
        end
    end

    if bD_max>abs(bD_min)
        if auc_nogasp(1)<0
            j=1;
        else
            j=2;
        end
    else
        if auc_nogasp(1)>0
            j=1;
        else
            j=2;
        end
    end
    while j+1<=length(auc_nogasp)
        TV_auc_nogasp(end+1)=abs(auc_nogasp(j))+abs(auc_nogasp(j+1));
        j=j+2;
    end
    TV_auc_nogasp_baseline=mean(TV_auc_nogasp);

    bD=breathData(baseline_ms(end)+1:end);
    TV_auc=[];
    auc=[];
    bD_max = max(bD);
    bD_min = min(bD);
    zc_indices = find(bD(1:end-1) .* bD(2:end) < 0);
      if length(sigh_idx)>0
        for i = 1:length(zc_indices)-1
            idx_start = zc_indices(i);
            idx_end = zc_indices(i+1);
            if idx_end<sigh_idx(checksigh)
                auc_nogasp(i) = sum(bD(idx_start:idx_end));
            elseif idx_start>sigh_idx(checksigh)
                auc_nogasp(i) = sum(bD(idx_start:idx_end));
            elseif (idx_end>sigh_idx(checksigh)) & (idx_start<sigh_idx(checksigh)) %exclude sigh
                if checksigh < length(sigh_idx)
                    checksigh = checksigh + 1;
                end
                continue
            end
            % Breath volume change
        end
    else
       
        for i = 1:length(zc_indices)-1
            idx_start = zc_indices(i);
             idx_end = zc_indices(i+1);
             auc_nogasp(i) = sum(breathData(idx_start:idx_end));
        end
    end
    
    
    if bD_max>abs(bD_min)
        if auc_nogasp(1)<0
            j=1;
        else
            j=2;
        end
    else
        if  auc_nogasp(1)>0
            j=1;
        else
            j=2;
        end
    end
    while j+1<=length(auc_nogasp)
        TV_auc_nogasp(end+1)=abs(auc_nogasp(j))+abs(auc_nogasp(j+1));
        j=j+2;
    end
    TV_auc_ci_nogasp=(TV_auc_nogasp-TV_auc_nogasp_baseline)./(TV_auc_nogasp+TV_auc_nogasp_baseline);

% % close all;
% 
% pbreath(pbreath>sighThd) = NaN;         % exclude sighs from TV analysis
% pbreathlow(pbreath>sighThd) = NaN;      % exclude sighs from TV analysis
% avepbreath = mean(pbreath,'omitnan');
% avepbreathlow = mean(pbreathlow,'omitnan');
% tidalVolume = avepbreath - avepbreathlow;
% breathRate = length(pbreath)/(length(timePeriod) / sampleRate);
% 

end

