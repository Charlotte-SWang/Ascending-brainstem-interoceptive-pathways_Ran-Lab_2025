function [breathRate,minuteVolume,tidalVolume,tidalVolume_auc,tidalVolume_auc_nogasp, minuteVolume_nogasp,nSigh, abspbthd,abspbthdlow, breathEvents] = analyzeBreath(ifplot,protocol_index, k, nn,breathData,timePeriod,sampleRate,maxBreathRate,pbreaththreshold,pbreaththresholdlow,pSighThd, abspbthreshold,abspbthresholdlow)
% Ran, analysis of breath rate and tidal volume
% Modified by Shiqi
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
junkpbreathlow = breathData(lk3); % Used to be lk3
avepbreath = mean(junkpbreath);
avepbreathlow = mean(junkpbreathlow);
tidalVolume = avepbreath - avepbreathlow;

    
%% analyzing breath rate and tidal volume
% if nargin == 9
    % abspbreaththreshold = abspbthreshold;
    % abspbreaththresholdlow = abspbthresholdlow;
% else
%     abspbreaththreshold = avepbreath - tidalVolume * pbreaththreshold;
%     abspbreaththresholdlow = avepbreathlow + tidalVolume * pbreaththresholdlow;
% end

breathEvents = [];
% sigh
%sighThd = avepbreathlow + tidalVolume * pSighThd;
sighThd = pSighThd;
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
if max(breathData) < abspbreaththreshold
    pbreath = 0;
    pbreathlow = 0;
else
    lk6 = [];
    [pbreath,breathEvents] = findpeaks(breathData,timePeriod,'MinPeakDistance',minPeakDis,'MinPeakHeight',abspbreaththreshold);
    
    invpresData = -breathData;
    [junk, lk5] = findpeaks(invpresData,timePeriod,'MinPeakDistance',minPeakDis,'MinPeakHeight',-abspbreaththresholdlow);
    for iInd = 1:length(lk5)   
        lk6(iInd) = find(timePeriod == lk5(iInd));
    end;
    if isempty(lk6) == 1
        pbreathlow = 0;
    else
        pbreathlow = breathData(lk6);
    if ifplot==1
%         plotting data
        figure;
        plot(timePeriod,breathData,breathEvents,pbreath,'ro');
%         plot(timePeriod,breathData);
%         line([timePeriod(1) timePeriod(end)],[abspbreaththresholdlow abspbreaththresholdlow],'color',[0.7 0.7 0.7],'linestyle','--');
%         line([timePeriod(1) timePeriod(end)],[abspbreaththreshold abspbreaththreshold],'color',[0.7 0.7 0.7],'linestyle','-.');
        hold on;
        plot(timePeriod,breathData,lk5,pbreathlow,'r*');
        title(['protocol' num2str(protocol_index), ' mouse',num2str(k),' trial',num2str(nn),' breathing analysis']);
    end
%         % use the following lines when plotting a part of the plot
%         timeLength = length(timePeriod);
%         timeLimit = floor(timeLength/10);
%         tempbreathEvents = breathEvents(breathEvents<timePeriod(timeLimit));
%         templk5 = lk5(lk5<timePeriod(timeLimit));
%         temppbreath = pbreath(1:length(tempbreathEvents));
%         temppbreathlow = pbreathlow(1:length(templk5));
%         figure;
%         plot(timePeriod(1:timeLimit),breathData(1:timeLimit),tempbreathEvents,temppbreath,'ro');
%         hold on;
%         plot(timePeriod(1:timeLimit),breathData(1:timeLimit),templk5,temppbreathlow,'r*');
%         title('breathing analysis');
    end;
end;
% close all;

sigh_idx=breathEvents(find(pbreath>sighThd));% exclude sighs from TV analysis 


closest_smaller = nan(size(sigh_idx)); % preallocate with NaN (if none found)

for k = 1:length(sigh_idx)
    smaller_vals = lk5(lk5 < sigh_idx(k));  % only smaller values
    if ~isempty(smaller_vals)
        closest_smaller(k) = max(smaller_vals); % pick the largest of the smaller ones
    end
end
sigh_lowidx = unique(closest_smaller);

pbreath(pbreath>sighThd) = NaN;         % exclude sighs from TV analysis
pbreathlow(pbreath>sighThd) = NaN;      % exclude sighs from TV analysis
avepbreath = mean(pbreath,'omitnan');
avepbreathlow = mean(pbreathlow,'omitnan');
tidalVolume = avepbreath - avepbreathlow;
breathRate = length(pbreath)/(length(timePeriod) / sampleRate);

TV_auc=[];
auc=[];
bD_max = max(breathData);
bD_min = min(breathData);
zc_indices = find(breathData(1:end-1) .* breathData(2:end) < 0);
sigh_idx = sort([sigh_idx sigh_lowidx]);
checksigh=1;

%%% Include gasp
for i = 1:length(zc_indices)-1
    idx_start = zc_indices(i);
    idx_end = zc_indices(i+1);
    auc(i) = sum(breathData(idx_start:idx_end));
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
tidalVolume_auc = mean(TV_auc);
minuteVolume =  sum(TV_auc);
checksigh=1
%%% Exclude gasp
if length(sigh_idx)>0
    for i = 1:length(zc_indices)-1
        idx_start = zc_indices(i);
        idx_end = zc_indices(i+1);
        if idx_end<sigh_idx(checksigh)
            auc_nogasp(i) = sum(breathData(idx_start:idx_end));
        elseif idx_start>sigh_idx(checksigh)
            auc_nogasp(i) = sum(breathData(idx_start:idx_end));
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

% if bD_max>abs(bD_min)
%     if auc(1)<0
%         j=1;
%     else
%         j=2;
%     end
% else
%     if auc(1)>0
%         j=1;
%     else
%         j=2;
%     end
% end
% while j+1<=length(auc)
%     TV_auc(end+1)=abs(auc(j))+abs(auc(j+1));
%     j=j+2;
% end

TV_auc_inhale=[];
TV_auc_exhale=[];
for j = 1:length(auc_nogasp)
    if bD_max*auc_nogasp(j)<0
      TV_auc_inhale(end+1)=abs(auc_nogasp(j))
    else 
      TV_auc_exhale(end+1)=abs(auc_nogasp(j));
    end
end
tidalVolume_auc_nogasp = mean(TV_auc_inhale)+mean(TV_auc_exhale);
minuteVolume_nogasp =  sum(TV_auc_inhale)+sum(TV_auc_exhale);
end