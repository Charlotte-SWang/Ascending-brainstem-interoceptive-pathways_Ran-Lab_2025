function [BP,SBP,DBP,PulsePressure] = analyzeBP(protocol_index, k, nn,BPData,timePeriod,sampleRate,EstMaxBP,MinPeakH_max,MinPeakH_min)
% Shiqi, analysis of Blood Pressure

    minPeakDis = 1/EstMaxBP*sampleRate;

    [sys,SBPevents] = findpeaks(BPData,timePeriod,'MinPeakDistance',minPeakDis,'MinPeakHeight',MinPeakH_max);
    invpresData = -BPData;
    [dia_minus, DBPevents] = findpeaks(invpresData,timePeriod,'MinPeakDistance',minPeakDis,'MinPeakHeight',-MinPeakH_min);
    dia = abs(dia_minus);
    if k==9
        figure;
        plot(timePeriod,BPData,SBPevents,sys,'ro');
        hold on;
        plot(timePeriod,BPData,DBPevents,dia,'r*');
        title(['protocol', num2str(protocol_index), ' mouse',num2str(k),' trial',num2str(nn),'BP analysis']);
    end
    SBP = mean(sys);
    DBP = mean(dia);
    BP = (SBP+DBP)/2;
    PulsePressure = SBP - DBP;
end

