function [nSigh] = analyzeSigh(ifplot,protocol_index, k, nn,breathData,timePeriod,sampleRate,maxSighRate,pSighThd,ifInvert)
 % Shiqi, analysis of Sigh
    minPeakDis = 1/maxSighRate*sampleRate;
    sighThd = pSighThd;
    if ifInvert == 1
        breathData = -breathData;
    end
    [psigh,lk4] = findpeaks(breathData,timePeriod,'MinPeakDistance',minPeakDis,'MinPeakHeight',sighThd);
    
    if ifplot==1
        figure;
        plot(timePeriod,breathData,lk4,psigh,'ko');
        line([timePeriod(1) timePeriod(end)],[sighThd sighThd],'color',[0 0 0],'linestyle',':');
        title(['protocol',num2str(protocol_index),'mouse',num2str(k),' trial',num2str(nn),' Sigh analysis']);
    end
    
    nSigh = length(lk4);

end

