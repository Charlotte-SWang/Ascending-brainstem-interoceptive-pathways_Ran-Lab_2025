% Code for saving individual trace--Breathing 5s
[fileDir] = uigetdir( 'Output');
cd(fileDir)
time = -4.999:0.001:10;
cp= data_extract{3,2};
ci = trialinfo_extract{3,2};
for i = 1:size(cp,2)
    %mm = (i-1)*3+1 
    %plot(data_extract{1,4}(:,i))
    %if ismember('VLM',trialinfo_extract{3,2}{i,3}) & ismember('E',trialinfo_extract{3,2}{i,4})
        figure('Visible','off');
        plot(time,cp(:,i))
        hold on
        lowbound = min(cp(:,i))-0.05;
        topbound = max(cp(:,i))+0.05;
        fill([0,0,5,5],[lowbound,topbound,topbound,lowbound],'blue','FaceAlpha',0.05,'linestyle','none');

            %fill([10,10,11,11],[lowbound,topbound,topbound,lowbound],'yellow','FaceAlpha',0.2,'linestyle','none');
            %hold on      
        xlim([-4.999,10])
        ylim([lowbound,topbound])
        xlabel('Time(s)');
        fname =[ci{i,4},'-',ci{i,3},'-',ci{i,2},'-T',num2str(ci{i,5}),'-idx',num2str(i),'-BreathingRawTrace-5s(1ms)'];
        title(fname, fontsize=12);
        saveas(gcf,strcat([fname,'.jpeg']))
        print(gcf,'-vector','-dsvg',[fname,'.svg']) % svg
        %print(gcf,[fname,'.tiff'],'-dtiff', '-r600') % svg
    %end
    %matname =['Protocol',num2str(protocol_index),'_Breathing_plot_data.mat']
    %save (matname,'trialinfo_extract','single_z_amp_tv_mouse','single_z_amp_tv_trial','single_z_amp_tv','single_ci_amp_tv','BR_zscore','BR_CI','protocol_index', 'BreathRate', 'amp_stage','Delta_amp','CI_amp','amp_stage_tv','Delta_amp_tv','CI_amp_tv','amp_stage_raw','bin')  
end

%% Code for saving individual trace--ECG 5s
[fileDir] = uigetdir( 'Output');
cd(fileDir)
time = -4.999:0.001:10;
cp= data_extract{3,1};
ci = trialinfo_extract{3,1};
for i = 1:size(cp,2)
    %mm = (i-1)*3+1 
    %plot(data_extract{1,4}(:,i))
    %if ismember('VLM',trialinfo_extract{3,2}{i,3}) & ismember('E',trialinfo_extract{3,2}{i,4})
        figure('Visible','off');
        plot(time,cp(:,i))
        hold on
        lowbound = min(cp(:,i))-0.05;
        topbound = max(cp(:,i))+0.05;
        fill([0,0,5,5],[lowbound,topbound,topbound,lowbound],'blue','FaceAlpha',0.05,'linestyle','none');

            %fill([10,10,11,11],[lowbound,topbound,topbound,lowbound],'yellow','FaceAlpha',0.2,'linestyle','none');
            %hold on      
        xlim([-4.999,10])
        ylim([lowbound,topbound])
        xlabel('Time(s)');
        fname =[ci{i,4},'-',ci{i,3},'-',ci{i,2},'-T',num2str(ci{i,5}),'-idx',num2str(i),'-ECGRawTrace-5s(1ms)'];
        title(fname, fontsize=12);
        saveas(gcf,strcat([fname,'.jpeg']))
        print(gcf,'-vector','-dsvg',[fname,'.svg']) % svg
        %print(gcf,[fname,'.tiff'],'-dtiff', '-r600') % svg
    %end
    %matname =['Protocol',num2str(protocol_index),'_Breathing_plot_data.mat']
    %save (matname,'trialinfo_extract','single_z_amp_tv_mouse','single_z_amp_tv_trial','single_z_amp_tv','single_ci_amp_tv','BR_zscore','BR_CI','protocol_index', 'BreathRate', 'amp_stage','Delta_amp','CI_amp','amp_stage_tv','Delta_amp_tv','CI_amp_tv','amp_stage_raw','bin')  
end
%% Code for saving individual trace--ECG 1min
[fileDir] = uigetdir( 'Output');
cd(fileDir)
xstart =-59.999;
xend = 63;
time = xstart:0.001:xend;
cp= data_extract{2,1};
ci = trialinfo_extract{2,1};
for i = 1:size(cp,2)
    %mm = (i-1)*3+1 
    %plot(data_extract{1,4}(:,i))
    %if ismember('VLM',trialinfo_extract{3,2}{i,3}) & ismember('E',trialinfo_extract{3,2}{i,4})
        figure('Visible','off');
        plot(time,cp(:,i))
        hold on
        lowbound = min(cp(:,i))-0.05;
        topbound = max(cp(:,i))+0.05;
        fill([0,0,60,60],[lowbound,topbound,topbound,lowbound],'blue','FaceAlpha',0.05,'linestyle','none');

            %fill([10,10,11,11],[lowbound,topbound,topbound,lowbound],'yellow','FaceAlpha',0.2,'linestyle','none');
            %hold on      
        xlim([xstart,xend])
        ylim([lowbound,topbound])
        xlabel('Time(s)');
        fname =[ci{i,4},'-',ci{i,3},'-',ci{i,2},'-T',num2str(ci{i,5}),'-idx',num2str(i),'-ECGRawTrace-1min(1ms)'];
        title(fname, fontsize=12);
        saveas(gcf,strcat([fname,'.jpeg']))
        print(gcf,'-vector','-dsvg',[fname,'.svg']) % svg
        %print(gcf,[fname,'.tiff'],'-dtiff', '-r600') % svg
    %end
    %matname =['Protocol',num2str(protocol_index),'_Breathing_plot_data.mat']
    %save (matname,'trialinfo_extract','single_z_amp_tv_mouse','single_z_amp_tv_trial','single_z_amp_tv','single_ci_amp_tv','BR_zscore','BR_CI','protocol_index', 'BreathRate', 'amp_stage','Delta_amp','CI_amp','amp_stage_tv','Delta_amp_tv','CI_amp_tv','amp_stage_raw','bin')  
end
%% Code for saving individual trace--BP 1min
[fileDir] = uigetdir( 'Output');
cd(fileDir)
xstart =-59.999;
xend = 63;
time = xstart:0.001:xend;
cp= data_extract{2,3};
ci = trialinfo_extract{2,3};
for i = 1:size(cp,2)
    %mm = (i-1)*3+1 
    %plot(data_extract{1,4}(:,i))
    %if ismember('VLM',trialinfo_extract{3,2}{i,3}) & ismember('E',trialinfo_extract{3,2}{i,4})
        figure('Visible','off');
        plot(time,cp(:,i))
        hold on
        lowbound = min(cp(:,i))-0.05;
        topbound = max(cp(:,i))+0.05;
        fill([0,0,60,60],[lowbound,topbound,topbound,lowbound],'blue','FaceAlpha',0.05,'linestyle','none');

            %fill([10,10,11,11],[lowbound,topbound,topbound,lowbound],'yellow','FaceAlpha',0.2,'linestyle','none');
            %hold on      
        xlim([xstart,xend])
        ylim([lowbound,topbound])
        xlabel('Time(s)');
        fname =[ci{i,4},'-',ci{i,3},'-',ci{i,2},'-T',num2str(ci{i,5}),'-idx',num2str(i),'-BPRawTrace-1min(1ms)'];
        title(fname, fontsize=12);
        saveas(gcf,strcat([fname,'.jpeg']))
        print(gcf,'-vector','-dsvg',[fname,'.svg']) % svg
        %print(gcf,[fname,'.tiff'],'-dtiff', '-r600') % svg
    %end
    %matname =['Protocol',num2str(protocol_index),'_Breathing_plot_data.mat']
    %save (matname,'trialinfo_extract','single_z_amp_tv_mouse','single_z_amp_tv_trial','single_z_amp_tv','single_ci_amp_tv','BR_zscore','BR_CI','protocol_index', 'BreathRate', 'amp_stage','Delta_amp','CI_amp','amp_stage_tv','Delta_amp_tv','CI_amp_tv','amp_stage_raw','bin')  
end
%% Code for saving individual trace--GP 5min
[fileDir] = uigetdir( 'Output');
cd(fileDir)
xstart =-299.999;
xend = 315;
time = xstart:0.001:xend;
stim = 300;
cp= data_extract{1,4};
ci = trialinfo_extract{1,4};
for i = 1:size(cp,2)
    %mm = (i-1)*3+1 
    %plot(data_extract{1,4}(:,i))
    %if ismember('VLM',trialinfo_extract{3,2}{i,3}) & ismember('E',trialinfo_extract{3,2}{i,4})
        figure('Visible','off');
        plot(time,cp(:,i))
        hold on
        lowbound = min(cp(:,i))-0.05;
        topbound = max(cp(:,i))+0.05;
        fill([0,0,stim,stim],[lowbound,topbound,topbound,lowbound],'blue','FaceAlpha',0.05,'linestyle','none');

            %fill([10,10,11,11],[lowbound,topbound,topbound,lowbound],'yellow','FaceAlpha',0.2,'linestyle','none');
            %hold on      
        xlim([xstart,xend])
        ylim([lowbound,topbound])
        xlabel('Time(s)');
        fname =[ci{i,4},'-',ci{i,3},'-',ci{i,2},'-T',num2str(ci{i,5}),'-idx',num2str(i),'-GPRawTrace-5min(1ms)'];
        title(fname, fontsize=12);
        saveas(gcf,strcat([fname,'.jpeg']))
        print(gcf,'-vector','-dsvg',[fname,'.svg']) % svg
        %print(gcf,[fname,'.tiff'],'-dtiff', '-r600') % svg
    %end
    %matname =['Protocol',num2str(protocol_index),'_Breathing_plot_data.mat']
    %save (matname,'trialinfo_extract','single_z_amp_tv_mouse','single_z_amp_tv_trial','single_z_amp_tv','single_ci_amp_tv','BR_zscore','BR_CI','protocol_index', 'BreathRate', 'amp_stage','Delta_amp','CI_amp','amp_stage_tv','Delta_amp_tv','CI_amp_tv','amp_stage_raw','bin')  
end

%% Code for saving individual trace--Breathing 5min
[fileDir] = uigetdir( 'Output');
cd(fileDir)
xstart =-299.999;
xend = 315;
time = xstart:0.001:xend;
stim = 300;
cp= data_extract{1,2};
ci = trialinfo_extract{1,2};
for i = 1:size(cp,2)
    %mm = (i-1)*3+1 
    %plot(data_extract{1,4}(:,i))
    %if ismember('VLM',trialinfo_extract{3,2}{i,3}) & ismember('E',trialinfo_extract{3,2}{i,4})
        %figure('Visible','off');
        figure;
        plot(time,cp(:,i))
        hold on
        lowbound = min(cp(:,i))-0.05;
        topbound = max(cp(:,i))+0.05;
        fill([0,0,stim,stim],[lowbound,topbound,topbound,lowbound],'blue','FaceAlpha',0.05,'linestyle','none');

            %fill([10,10,11,11],[lowbound,topbound,topbound,lowbound],'yellow','FaceAlpha',0.2,'linestyle','none');
            %hold on      
        xlim([xstart,xend])
        ylim([lowbound,topbound])
        xlabel('Time(s)');
        fname =[ci{i,4},'-',ci{i,3},'-',ci{i,2},'-T',num2str(ci{i,5}),'-idx',num2str(i),'-BreathingRawTrace-5min(1ms)'];
        title(fname, fontsize=12);
        saveas(gcf,strcat([fname,'.jpeg']))
        print(gcf,'-vector','-dsvg',[fname,'.svg']) % svg
        %print(gcf,[fname,'.tiff'],'-dtiff', '-r600') % svg
    %end
    %matname =['Protocol',num2str(protocol_index),'_Breathing_plot_data.mat']
    %save (matname,'trialinfo_extract','single_z_amp_tv_mouse','single_z_amp_tv_trial','single_z_amp_tv','single_ci_amp_tv','BR_zscore','BR_CI','protocol_index', 'BreathRate', 'amp_stage','Delta_amp','CI_amp','amp_stage_tv','Delta_amp_tv','CI_amp_tv','amp_stage_raw','bin')  
end

%% %%%%%%%%%% Close all figures 
allFigs = findall(0, 'Type', 'figure');
allFigs(allFigs == gcf) = [];
close(allFigs);

%% CLOSE INVISIBLE FIGURES
figs = findall(0, 'Type', 'figure', 'Visible', 'off');
close(figs);


%% ONLY CLOSE INVISIBLE FIGURES BUT KEEP VISIBLE OPENS
allFigs = findall(0, 'Type', 'figure');
for i = 1:numel(allFigs)
    if strcmp(allFigs(i).Visible, 'off')
        close(allFigs(i));
    end
end



%% Plotting individual trace--Breathing 5min--detrend 

%%%   PLOT ONLY !!!!!!
[fileDir] = uigetdir( 'Output');
cd(fileDir)
% xstart =-299.999;
% xend = 315;
% time = xstart:0.001:xend;
% stim = 300;

xstart =-4.999;
xend = 10;
time = xstart:0.001:xend;
stim = 5;
cp= data_extract{3,2};
ci = trialinfo_extract{3,2};
for i = 1:size(cp,2)
    %mm = (i-1)*3+1 
    %plot(data_extract{1,4}(:,i))
    %if ismember('VLM',trialinfo_extract{3,2}{i,3}) & ismember('E',trialinfo_extract{3,2}{i,4})
        %figure('Visible','off');
        figure;
        plot(time,detrend(cp(:,i)))
        hold on
        lowbound = min(detrend(cp(:,i)))-0.05;
        topbound = max(detrend(cp(:,i)))+0.05;
        fill([0,0,stim,stim],[lowbound,topbound,topbound,lowbound],'blue','FaceAlpha',0.05,'linestyle','none');

            %fill([10,10,11,11],[lowbound,topbound,topbound,lowbound],'yellow','FaceAlpha',0.2,'linestyle','none');
            %hold on      
        xlim([xstart,xend])
        ylim([lowbound,topbound])
        xlabel('Time(s)');
        fname =[ci{i,4},'-',ci{i,3},'-',ci{i,2},'-T',num2str(ci{i,5}),'-idx',num2str(i),'-BreathingRawTrace-5min(1ms)'];
        title(fname, fontsize=12);
        %saveas(gcf,strcat([fname,'.jpeg']))
        %print(gcf,'-vector','-dsvg',[fname,'.svg']) % svg
        %print(gcf,[fname,'.tiff'],'-dtiff', '-r600') % svg
    %end
    %matname =['Protocol',num2str(protocol_index),'_Breathing_plot_data.mat']
    %save (matname,'trialinfo_extract','single_z_amp_tv_mouse','single_z_amp_tv_trial','single_z_amp_tv','single_ci_amp_tv','BR_zscore','BR_CI','protocol_index', 'BreathRate', 'amp_stage','Delta_amp','CI_amp','amp_stage_tv','Delta_amp_tv','CI_amp_tv','amp_stage_raw','bin')  
end