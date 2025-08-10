 %% plot GP
j=1;
o=1;
l=1;
m=1;
jj=1;
oo=1;
ll=1;
mm=1;
TTV_a = [];
TTV_m = [];
TTV_p = [];
TTV_c = [];
mTTV_a = [];
mTTV_m = [];
mTTV_p = [];
mTTV_c = [];
cc=GP;
%cc =CI_amp_tv;
%cc = MAP_zscore
for k = 1:size(info,1)
    kk = (k-1)*nrep +1;
    if ismember('C',info{k,2}(:))
        mTTV_c(jj,:) = (cc{kk}+cc{kk+1}+cc{kk+2})/3;
        jj=jj+1;
    else
        if ismember('A',info(k,2))
            mTTV_a(oo,:) =   (cc{kk}+cc{kk+1}+cc{kk+2})/3;
            oo=oo+1;
        elseif ismember('M',info(k,2)) 
            mTTV_m(ll,:) =  (cc{kk}+cc{kk+1}+cc{kk+2})/3;
            ll=ll+1;
        elseif ismember('P',info(k,2))
            mTTV_p(mm,:) =   (cc{kk}+cc{kk+1}+cc{kk+2})/3;
            mm=mm+1;
        end
    end
    for nn = 1:nrep
        kk = (k-1)*nrep +nn;
        if ismember('C',info{k,2}(:))
            TTV_c(j,:) = cc{kk};
            j = j+1;
        else
            if ismember('A',info(k,2))
                TTV_a(o,:) = cc{kk};
                o = o+1;
            elseif ismember('M',info(k,2)) 
                TTV_m(l,:) = cc{kk};
                l = l+1;
            elseif ismember('P',info(k,2))
                TTV_p(m,:) = cc{kk};
                m = m+1;
            end
        end
    end
end 

% Heatmap plot total volumne_averaged
 % tt = 135;
% tt=80;
% bintime = 1;
% stim = 60; %stim=60;
% nbin = stim/bintime
% %st =(-nbin+0.5)*bintime;
% st =(-20+0.5)*bintime;
% et =tt-(20+0.5)*bintime;
% time = st:bintime:et
% figure;
% %clims =[-0.1,0.1];
% clims = [-1,1]
% showperiod = [41:120]

%10S--5S
tt=240
stim = 600000; %stim=60;
bintime = 1/sampleRate;
nbin = stim/bintime;
%st =(-nbin+0.5)*bintime;
st =(-60+0.5);
et =tt-(-60+0.5);
time = st:bintime:et
figure;
%clims =[-0.1,0.1];
clims = [8,15]
showperiod = [61*sampleRate:240*sampleRate]

colormap('jet')
tcl=tiledlayout(4,1);
nexttile
imagesc(mTTV_a([1:5 7:8],showperiod ),'XData',time,clims);
title('anterior');
nexttile
imagesc(mTTV_m([1:5 7:10],showperiod ),'XData',time,clims);
title('intermediate');
nexttile
imagesc(mTTV_p([1:7 9:13 15:19],showperiod ),'XData',time,clims);
title('posterior');
nexttile
imagesc(mTTV_c(1:end,showperiod),'XData',time,clims);
title('control');

cb = colorbar;
cb.Layout.Tile = 'east';
title(tcl,' GP mouse average 3 trials, bin=1s ', fontsize=12);
%% plot BREATHING RATE
nrep=3
j=1;
o=1;
l=1;
m=1;
jj=1;
oo=1;
ll=1;
mm=1;
TTV_a = [];
TTV_m = [];
TTV_p = [];
TTV_c = [];
mTTV_a = [];
mTTV_m = [];
mTTV_p = [];
mTTV_c = [];
%cc=BreathRate;
cc=BR_CI
%cc =CI_amp_tv;
%cc=amp_stage_raw;
for k = 1:size(info,1)
    kk = (k-1)*nrep +1;
    if ismember('C',info{k,2}(:))
        mTTV_c(jj,:) = (cc{kk}+cc{kk+1}+cc{kk+2})/3;
        jj=jj+1;
    else
        if ismember('A',info(k,2))
            mTTV_a(oo,:) =   (cc{kk}+cc{kk+1}+cc{kk+2})/3;
            oo=oo+1;
        elseif ismember('M',info(k,2)) 
            mTTV_m(ll,:) =  (cc{kk}+cc{kk+1}+cc{kk+2})/3;
            ll=ll+1;
        elseif ismember('P',info(k,2))
            mTTV_p(mm,:) =   (cc{kk}+cc{kk+1}+cc{kk+2})/3;
            mm=mm+1;
        end
    end
    for nn = 1:nrep
        kk = (k-1)*nrep +nn;
        if ismember('C',info{k,2}(:))
            TTV_c(j,:) = cc{kk};
            j = j+1;
        else
            if ismember('A',info(k,2))
                TTV_a(o,:) = cc{kk};
                o = o+1;
            elseif ismember('M',info(k,2)) 
                TTV_m(l,:) = cc{kk};
                l = l+1;
            elseif ismember('P',info(k,2))
                TTV_p(m,:) = cc{kk};
                m = m+1;
            end
        end
    end
end 

% Heatmap plot total volumne_averaged
 % tt = 135;
% tt=80;
% bintime = 1;
% stim = 60; %stim=60;
% nbin = stim/bintime
% %st =(-nbin+0.5)*bintime;
% st =(-20+0.5)*bintime;
% et =tt-(20+0.5)*bintime;
% time = st:bintime:et
% figure;
% %clims =[-0.1,0.1];
% clims = [-1,1]
% showperiod = [41:120]

%10S--5S
tt=20;
bintime = 1;
stim = 10; %stim=60;
nbin = stim/bintime
%st =(-nbin+0.5)*bintime;
st =(-10+0.5)*bintime;
et =tt-(10+0.5)*bintime;
time = st:bintime:et
figure;
clims =[-1 0];
%clims = [0,400]
showperiod = [1:20]

colormap('jet')
tcl=tiledlayout(4,1);
nexttile
imagesc(mTTV_a([1:5 7:8],showperiod ),'XData',time,clims);
title('anterior');
nexttile
imagesc(mTTV_m([1:5 7:10],showperiod ),'XData',time,clims);
title('intermediate');
nexttile
imagesc(mTTV_p([1:7 9:13 15:19],showperiod ),'XData',time,clims);
title('posterior');
nexttile
imagesc(mTTV_c(1:end,showperiod),'XData',time,clims);
title('control');

cb = colorbar;
cb.Layout.Tile = 'east';
title(tcl,' 10S Breath Rate CI mouse average 3 trials, bin=1s ', fontsize=12);

%% Plot Breath Rate (bin =1s=)
    tt = 15;
    bintime = 1;
    stim = 10;
    nbin = stim/bintime;
    st =(-5)*bintime;
    et =tt-(5)*bintime;
    time = st:bintime:et;
    showperiod = [5:20]
    % ylim([lowbound,topbound])
    % xlim([(-nbin+1)*bintime (tt/bin-nbin)*bintime])
    % set(gca,'FontSize', 15);
    % hold on
    % 
    % fill([0,0,120,120],[lowbound,topbound,topbound,lowbound],'yellow','FaceAlpha',0.2,'linestyle','none');
    % 
    lowbound = 0
    topbound = 400
    exchange_list='0123456789ABCDEF#'
    %colors ={'#0072BD','#D95319','#EDB120',	'#7E2F8E',	'#77AC30','#4DBEEE','#A2142F','#FF69B4'}
    colors = { '#E274A9', '#13AF68','#0F80FF','#424242'};  % line colors %"#E274A9" "#26A7E1'#FFE009'

    x=time
    h = figure('visible','on');
    hA = axes(h);

    mm=[]
    error=[]
    var = mTTV_a([1:5 7:8],showperiod );
    mm(:,1) = mean(var,'omitnan')
    error(:,1)= std(var,'omitnan')/sqrt(size(var,1));
    var = mTTV_m([1:5 7:10],showperiod );
    mm(:,2) = mean(var,'omitnan')
    error(:,2)= std(var,'omitnan')/sqrt(size(var,1));
    var = mTTV_p([1:7 9:13 15:19],showperiod );
    mm(:,3) = mean(var,'omitnan')
    error(:,3)= std(var,'omitnan')/sqrt(size(var,1));
    var = mTTV_c(:,showperiod );
    mm(:,4) = mean(var,'omitnan')
    error(:,4)= std(var,'omitnan')/sqrt(size(var,1));
    count=1
    p=[]
    for i=1:size(mm,2)
        y1=mm(:,i)+error(:,i)
        y2=mm(:,i)-error(:,i)
        x2=[x,fliplr(x)];
        inBetween = [y1',fliplr(y2')];
        string=(colors{count})
        num=zeros(1,3)
        for r=1:3
            temp_Coe1=find(exchange_list==string(r*2))-1;
            temp_Coe2=find(exchange_list==string(r*2+1))-1;
            num(r)=(16*temp_Coe1+temp_Coe2)/255;
        end
        fi=fill(x2,inBetween,num,'HandleVisibility','off');
        set(fi,'edgealpha',0,'facealpha',0.3)
        hold on
        plot(x,mm(:,i), 'LineWidth',1, 'Color', colors{count});

        count=count+1
        ylim([lowbound,topbound])
        xlim([st et])
        set(gca,'FontSize', 15);
    end
    hold on 
    ylims = ylim;  % Get current y-axis limits
    line([0 0], ylims, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1.5);
    xlims = xlim;  % Get current y-axis limits
    hold on
    %line(xlims,[0 0],  'Color', 'k', 'LineStyle', '--', 'LineWidth', 1.5);
   
    group = {'Anterior','intermediate','posterior','Control'}
    lgd=legend(group,'Location','NorthEastOutside','Box','off')
    title(lgd,'Group')
    
    %fill([0,0,stim,stim],[lowbound,topbound,topbound,lowbound],'yellow','FaceAlpha',0.2,'linestyle','none');


%% plot BREATHING RATE_change index
nrep=3
j=1;
o=1;
l=1;
m=1;
jj=1;
oo=1;
ll=1;
mm=1;
TTV_a = [];
TTV_m = [];
TTV_p = [];
TTV_c = [];
mTTV_a = [];
mTTV_m = [];
mTTV_p = [];
mTTV_c = [];
cc=BR_CI;
%cc =CI_amp_tv;
%cc=amp_stage_raw;
for k = 1:size(info,1)
    kk = (k-1)*nrep +1;
    if ismember('C',info{k,2}(:))
        mTTV_c(jj,:) = (cc{kk}+cc{kk+1}+cc{kk+2})/3;
        jj=jj+1;
    else
        if ismember('A',info(k,2))
            mTTV_a(oo,:) =   (cc{kk}+cc{kk+1}+cc{kk+2})/3;
            oo=oo+1;
        elseif ismember('M',info(k,2)) 
            mTTV_m(ll,:) =  (cc{kk}+cc{kk+1}+cc{kk+2})/3;
            ll=ll+1;
        elseif ismember('P',info(k,2))
            mTTV_p(mm,:) =   (cc{kk}+cc{kk+1}+cc{kk+2})/3;
            mm=mm+1;
        end
    end
    for nn = 1:nrep
        kk = (k-1)*nrep +nn;
        if ismember('C',info{k,2}(:))
            TTV_c(j,:) = cc{kk};
            j = j+1;
        else
            if ismember('A',info(k,2))
                TTV_a(o,:) = cc{kk};
                o = o+1;
            elseif ismember('M',info(k,2)) 
                TTV_m(l,:) = cc{kk};
                l = l+1;
            elseif ismember('P',info(k,2))
                TTV_p(m,:) = cc{kk};
                m = m+1;
            end
        end
    end
end 

% Heatmap plot total volumne_averaged
 % tt = 135;
% tt=80;
% bintime = 1;
% stim = 60; %stim=60;
% nbin = stim/bintime
% %st =(-nbin+0.5)*bintime;
% st =(-20+0.5)*bintime;
% et =tt-(20+0.5)*bintime;
% time = st:bintime:et
% figure;
 %clims =[-0.1,0.1];
% clims = [-1,1]
% showperiod = [41:120]

%10S--5S
tt=20;
bintime = 1;
stim = 10; %stim=60;
nbin = stim/bintime
%st =(-nbin+0.5)*bintime;
st =(-10+0.5)*bintime;
et =tt-(10+0.5)*bintime;
time = st:bintime:et
figure;
%clims =[-0.2,0.2];
clims = [-1,1]
showperiod = [1:20]

colormap('jet')
tcl=tiledlayout(4,1);
nexttile
imagesc(mTTV_a([1:5 7:8],showperiod ),'XData',time,clims);
title('anterior');
nexttile
imagesc(mTTV_m([1:5 7:10],showperiod ),'XData',time,clims);
title('intermediate');
nexttile
imagesc(mTTV_p([1:7 9:13 15:19],showperiod ),'XData',time,clims);
title('posterior');
nexttile
imagesc(mTTV_c(1:end,showperiod),'XData',time,clims);
title('control');

cb = colorbar;
cb.Layout.Tile = 'east';
title(tcl,' 10S TV CI mouse average 3 trials, bin=1s ', fontsize=12);

    %% Plot Breath Rate Normalization (bin =1s=)
    tt = 15;
    bintime = 1;
    stim = 10;
    nbin = stim/bintime;
    st =(-5)*bintime;
    et =tt-(5)*bintime;
    time = st:bintime:et;
    showperiod = [5:20]
    % ylim([lowbound,topbound])
    % xlim([(-nbin+1)*bintime (tt/bin-nbin)*bintime])
    % set(gca,'FontSize', 15);
    % hold on
    % 
    % fill([0,0,120,120],[lowbound,topbound,topbound,lowbound],'yellow','FaceAlpha',0.2,'linestyle','none');
    % 
    lowbound = -0.6
    topbound = 0.15
    exchange_list='0123456789ABCDEF#'
    %colors ={'#0072BD','#D95319','#EDB120',	'#7E2F8E',	'#77AC30','#4DBEEE','#A2142F','#FF69B4'}
    colors = { '#E274A9', '#13AF68','#0F80FF','#424242'};  % line colors %"#E274A9" "#26A7E1'#FFE009'

    x=time
    h = figure('visible','on');
    hA = axes(h);

    mm=[]
    error=[]
    var = mTTV_a([1:5 7:8],showperiod );
    mm(:,1) = mean(var,'omitnan')
    error(:,1)= std(var,'omitnan')/sqrt(size(var,1));
    var = mTTV_m([1:5 7:10],showperiod );
    mm(:,2) = mean(var,'omitnan')
    error(:,2)= std(var,'omitnan')/sqrt(size(var,1));
    var = mTTV_p([1:7 9:13 15:19],showperiod );
    mm(:,3) = mean(var,'omitnan')
    error(:,3)= std(var,'omitnan')/sqrt(size(var,1));
    var = mTTV_c(:,showperiod );
    mm(:,4) = mean(var,'omitnan')
    error(:,4)= std(var,'omitnan')/sqrt(size(var,1));
    count=1
    p=[]
    for i=1:size(mm,2)
        y1=mm(:,i)+error(:,i)
        y2=mm(:,i)-error(:,i)
        x2=[x,fliplr(x)];
        inBetween = [y1',fliplr(y2')];
        string=(colors{count})
        num=zeros(1,3)
        for r=1:3
            temp_Coe1=find(exchange_list==string(r*2))-1;
            temp_Coe2=find(exchange_list==string(r*2+1))-1;
            num(r)=(16*temp_Coe1+temp_Coe2)/255;
        end
        fi=fill(x2,inBetween,num,'HandleVisibility','off');
        set(fi,'edgealpha',0,'facealpha',0.3)
        hold on
        plot(x,mm(:,i), 'LineWidth',1, 'Color', colors{count});

        count=count+1
        ylim([lowbound,topbound])
        xlim([st et])
        set(gca,'FontSize', 15);
    end
    hold on 
    ylims = ylim;  % Get current y-axis limits
    line([0 0], ylims, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1.5);
    xlims = xlim;  % Get current y-axis limits
    hold on
    line(xlims,[0 0],  'Color', 'k', 'LineStyle', '--', 'LineWidth', 1.5);
   
    group = {'Anterior','intermediate','posterior','Control'}
    lgd=legend(group,'Location','NorthEastOutside','Box','off')
    title(lgd,'Group')
    
    %fill([0,0,stim,stim],[lowbound,topbound,topbound,lowbound],'yellow','FaceAlpha',0.2,'linestyle','none');


%% plot BREATHING tidal volume--preparation
nrep=3
j=1;
o=1;
l=1;
m=1;
jj=1;
oo=1;
ll=1;
mm=1;
TTV_a = [];
TTV_m = [];
TTV_p = [];
TTV_c = [];
mTTV_a = {};
mTTV_m = {};
mTTV_p = {};
mTTV_c = {};

gmTTV{1}=[];
gmTTV{2}=[];
gmTTV{3}=[];
gmTTV{4}=[];
cc =single_ci_amp_tv;
groupinfo=[];
imTTV={};
for k = 1:size(info,1)
    kk = (k-1)*nrep +1;
    if ismember('C',info{k,2}(:)) 
        mTTV_c{jj} = [cc{kk}; cc{kk+1}; cc{kk+2}];
        gmTTV{4}=[gmTTV{4}; mTTV_c{jj}];
       
        imTTV=[imTTV; mTTV_c{jj}];
        groupinfo=[groupinfo 4];

        jj=jj+1;
        j=j+1;
    else
        if ismember('A',info(k,2)) 
            if o ~=6
                mTTV_a{oo} = [cc{kk}; cc{kk+1}; cc{kk+2}];
                gmTTV{1}=[gmTTV{1}; mTTV_a{oo}];
    
                imTTV=[imTTV; mTTV_a{oo}];
    
                groupinfo=[groupinfo 1];
    
                oo=oo+1;
            end
            o=o+1;
        elseif ismember('M',info(k,2))
            if l ~=6
                mTTV_m{ll} = [cc{kk}; cc{kk+1}; cc{kk+2}];
                gmTTV{2}=[gmTTV{2}; mTTV_m{ll}];
    
                imTTV=[imTTV; mTTV_m{ll}];
                groupinfo=[groupinfo 2];
    
                ll=ll+1;
            end
            l=l+1;
        elseif ismember('P',info(k,2))
            if  (m ~=8) & (m ~=14)
                mTTV_p{mm} = [cc{kk}; cc{kk+1}; cc{kk+2}];
                gmTTV{3}=[gmTTV{3}; mTTV_p{mm}];

                imTTV=[imTTV; mTTV_p{mm}];
                groupinfo=[groupinfo 3];

                mm=mm+1;
            end
            m=m+1;
        end
    end
end

%% plot BREATHING tidal volume-group 3d histogram
figure
xs = 0:5:15
colors = jet(numel(xs));
for i = 1:numel(xs)
    histogram2(repmat(xs(i),1,size(gmTTV{i},1)),gmTTV{i}.','FaceColor',colors(i,:),'Normalization','percentage')
    hold on
end
hold off

%% plot BREATHING tidal volume- individual 3d histogram
figure
xs = 1:1:44
colormap('jet')
colors = jet(numel(groupinfo));
for i = 1:numel(xs)
%    histogram2(repmat(xs(i),1,size(imTTV{i},1)),imTTV{i}.','BinWidth',[1,0.02],'FaceColor',colors(i,:),'Normalization','percentage')
    histogram2(repmat(xs(i),1,size(imTTV{i},1)),imTTV{i}.','BinWidth',[0.05,0.01],'FaceColor','flat','Normalization','percentage')
    hold on
end
colorbar
%caxis([0 30])
hold off

%% plot BREATHING tidal volume- individual 2d heatmap
figure
xs = 0:1:43
colormap('jet')

[groupinfo_sorted, rowOrder] = sort(groupinfo, 'ascend');
imTTV_sorted = imTTV(rowOrder);

for i = 1:4
    % Rank rows by their mean values (descending)
    idx = find(groupinfo_sorted == i)
    data = imTTV_sorted(idx,:)
    rowMeans2 = cellfun(@mean, data)
    [~, rowOrder2] = sort(rowMeans2, 'descend');
    imTTV_sorted2(idx,:) = data(rowOrder2, :);
    imTTV_group{i}=data(rowOrder2, :);
end

% Plot histogram2 heatmap
ax1 = axes('Position', [0.1, 0.25, 0.6, 0.65]);
axes(ax1); hold on
colormap(ax1, 'jet');
groupinfo_sorted=groupinfo(rowOrder)
for i = 1:numel(xs)
    histogram2(repmat(xs(i),1,size(imTTV_sorted2{i},1)),imTTV_sorted2{i}.','BinWidth',[1,0.05],'FaceColor','flat','Normalization','percentage','DisplayStyle','tile','ShowEmptyBins','on')
    hold on
end

hold off
xlim([min(xs), max(xs)+1]);  % ensure full axis range

caxis([0 20])


% Create a group indicator bar
%subplot(2,1,2);  % Right subplot for group bar
ax2 = axes('Position', [0.1, 0.1, 0.6, 0.05]);
imagesc(groupinfo_sorted);  % must be a row vector for vertical bar
colormap(gca, lines(numel(unique(groupinfo_sorted))));
set(gca, 'XTick', [], 'YTick', 1:length(groupinfo_sorted), 'YDir', 'normal');
ylabel('Groups');

ax3 = axes('Position', [0.7, 0.25, 0.1, 0.65]);  % slim colorbar to the right
%caxis(ax2, [min(imTTV_sorted(:)) max(data_sorted(:))]);  % match color limits
caxis(ax3,[0 20]);
colorbar(ax3)%, 'Location', 'eastoutside');  % Attach colorbar
axis off;  % hide dummy axis
%% plot BREATHING tidal volume- 2d histogram
figure

%colormap('jet')
subplot(4, 1, 1)
histogram(gmTTV{1}.','BinWidth',0.01,'Normalization','count','FaceColor','Red')
xlim([-0.2 0.5])
ylim([0 30])
subplot(4, 1, 2)
histogram(gmTTV{2}.','BinWidth',0.01,'Normalization','count','FaceColor','Green')
xlim([-0.2 0.5])
ylim([0 30])
subplot(4, 1, 3)
histogram(gmTTV{3}.','BinWidth',0.01,'Normalization','count','FaceColor','Blue')
xlim([-0.2 0.5])
ylim([0 30])
subplot(4, 1, 4)
histogram(gmTTV{4}.','BinWidth',0.01,'Normalization','count','FaceColor','Yellow')
xlim([-0.2 0.5])
ylim([0 20])
%caxis([0 100])
hold off

%% plot BREATHING tidal volume- Breathing Tidal Volume-Preview of Kernel Fitting
figure 
% Plot normalized histogram (as a density, not counts)
normaldata = gmTTV{2};
normalpd = fitdist(normaldata,"Kernel");
plot(normalpd);
%% plot BREATHING tidal volume- Breathing Tidal Volume-Final Cumulative distribution
%colors = { "#E274A9", '#13AF68',"#26A7E1","#E95412"};  % line colors %"#E274A9" "#26A7E1'#FFE009'
colors = { "#E274A9", '#13AF68',"#0F80FF","#424242"};  % line colors %"#E274A9" "#26A7E1'#FFE009'

hold on;

% Define common x-range
%x = linspace(min(cellfun(@min, gmTTV)), max(cellfun(@max, gmTTV)), 200);
x = linspace(-1,1,200)
for i = 1:4
    data = gmTTV{i};
    pd = fitdist(data, 'Kernel');
    y = cdf(pd, x);  % evaluate CDF
    plot(x, y, 'Color', colors{i}, 'LineWidth',2);  % clean line only
end

xlabel('Change index');
ylabel('Cumulative distribution');
legend('Anterior', 'Intermediate', 'Posterior', 'Control');
title('Overlaid CDFs from Kernel Distributions of 10s tidal volume data');
grid off;
%% plot BREATHING tidal volume-Breathing Tidal Volume-Density Cumulative distribution
figure
% Plot kernel density estimate
colors = { "#E274A9", '#13AF68',"#0F80FF","#424242"};  % line colors %"#E274A9" "#26A7E1'#FFE009'
%colors = { "#E95412", '#13AF68',"#26A7E1",'#FFE009'};  % line colors %"#E274A9" "#26A7E1
hold on;

% Define common x-range
%x = linspace(min(cellfun(@min, gmTTV)), max(cellfun(@max, gmTTV)), 200);
x = linspace(-1,0.6,200)
for i = 1:4
    [f, x] = ksdensity(gmTTV{i});
    plot(x, f, 'r-', 'LineWidth', 2,'Color',colors{i});
    hold on;
end

% Labels and formatting
xlabel('Change index');
ylabel('Density');
xlim([-1,0.6])
legend('Anterior', 'Intermediate', 'Posterior', 'Control');
title('Histogram with Kernel Density Estimate of 10s tidal volume data');


 %% plot Heart RATE
j=1;
o=1;
l=1;
m=1;
jj=1;
oo=1;
ll=1;
mm=1;
TTV_a = [];
TTV_m = [];
TTV_p = [];
TTV_c = [];
mTTV_a = [];
mTTV_m = [];
mTTV_p = [];
mTTV_c = [];
%cc=BR_CI;
cc = HeartRate;
%cc = MAP_zscore
for k = 1:size(info,1)
    kk = (k-1)*nrep +1;
    if ismember('C',info{k,2}(:))
        mTTV_c(jj,:) = (cc{kk}+cc{kk+1}+cc{kk+2})/3;
        jj=jj+1;
    else
        if ismember('A',info(k,2))
            mTTV_a(oo,:) =   (cc{kk}+cc{kk+1}+cc{kk+2})/3;
            oo=oo+1;
        elseif ismember('M',info(k,2)) 
            mTTV_m(ll,:) =  (cc{kk}+cc{kk+1}+cc{kk+2})/3;
            ll=ll+1;
        elseif ismember('P',info(k,2))
            mTTV_p(mm,:) =   (cc{kk}+cc{kk+1}+cc{kk+2})/3;
            mm=mm+1;
        end
    end
    for nn = 1:nrep
        kk = (k-1)*nrep +nn;
        if ismember('C',info{k,2}(:))
            TTV_c(j,:) = cc{kk};
            j = j+1;
        else
            if ismember('A',info(k,2))
                TTV_a(o,:) = cc{kk};
                o = o+1;
            elseif ismember('M',info(k,2)) 
                TTV_m(l,:) = cc{kk};
                l = l+1;
            elseif ismember('P',info(k,2))
                TTV_p(m,:) = cc{kk};
                m = m+1;
            end
        end
    end
end 

% Heatmap plot total volumne_averaged
 % tt = 135;
tt=1500;
bintime = 3;
stim = 120; %stim=60;
nbin = stim/bintime
%st =(-nbin+0.5)*bintime;
st =(-30+0.5)/bintime;
et =(tt-(30+0.5))/bintime;
time = st:bintime:et
figure;
clims =[400,800];
%clims = [-1,1]
showperiod = [31:150]

colormap('jet')
tcl=tiledlayout(4,1);
nexttile
imagesc(mTTV_a([1:5 7:8],showperiod ),'XData',time,clims);
title('anterior');
nexttile
imagesc(mTTV_m(1:end,showperiod ),'XData',time,clims);
title('intermediate');
nexttile
imagesc(mTTV_p(1:end,showperiod ),'XData',time,clims);
title('posterior');
nexttile
imagesc(mTTV_c(1:end,showperiod),'XData',time,clims);
title('control');

cb = colorbar;
cb.Layout.Tile = 'east';
title(tcl,' Heart rate change index mouse average 3 trials, bin=1s, stimulation 0-25s ', fontsize=12);

%% plot total volumne

j=1;
o=1;
l=1;
m=1;
jj=1;
oo=1;
ll=1;
mm=1;
TTV_a = [];
TTV_m = [];
TTV_p = [];
TTV_c = [];
mTTV_a = [];
mTTV_m = [];
mTTV_p = [];
mTTV_c = [];
cc=Delta_BP
%cc = MAP_zscore
for k = 1:size(info,1)
    kk = (k-1)*nrep +1;
    if ismember('C',info{k,2}(:))
        mTTV_c(jj,:) = (cc{kk}+cc{kk+1}+cc{kk+2})/3;
        jj=jj+1;
    else
        if ismember('A',info(k,2))
            mTTV_a(oo,:) =   (cc{kk}+cc{kk+1}+cc{kk+2})/3;
            oo=oo+1;
        elseif ismember('M',info(k,2)) 
            mTTV_m(ll,:) =  (cc{kk}+cc{kk+1}+cc{kk+2})/3;
            ll=ll+1;
        elseif ismember('P',info(k,2))
            mTTV_p(mm,:) =   (cc{kk}+cc{kk+1}+cc{kk+2})/3;
            mm=mm+1;
        end
    end
    for nn = 1:nrep
        kk = (k-1)*nrep +nn;
        if ismember('C',info{k,2}(:))
            TTV_c(j,:) = cc{kk};
            j = j+1;
        else
            if ismember('A',info(k,2))
                TTV_a(o,:) = cc{kk};
                o = o+1;
            elseif ismember('M',info(k,2)) 
                TTV_m(l,:) = cc{kk};
                l = l+1;
            elseif ismember('P',info(k,2))
                TTV_p(m,:) = cc{kk};
                m = m+1;
            end
        end
    end
 end

 %% plot BP value
nrep=3
j=1;
o=1;
l=1;
m=1;
jj=1;
oo=1;
ll=1;
mm=1;
TTV_a = [];
TTV_m = [];
TTV_p = [];
TTV_c = [];
mTTV_a = [];
mTTV_m = [];
mTTV_p = [];
mTTV_c = [];
cc = Delta_BP
for k = 1:size(info,1)
    kk = (k-1)*nrep +1;
    if ismember('C',info{k,2}(:))
        mTTV_c(jj,:) = (cc{kk}+cc{kk+1}+cc{kk+2})/3;
        jj=jj+1;
    else
        if ismember('A',info(k,2))
            mTTV_a(oo,:) =   (cc{kk}+cc{kk+1}+cc{kk+2})/3;
            oo=oo+1;
        elseif ismember('M',info(k,2)) 
            mTTV_m(ll,:) =  (cc{kk}+cc{kk+1}+cc{kk+2})/3;
            ll=ll+1;
        elseif ismember('P',info(k,2))
            mTTV_p(mm,:) =   (cc{kk}+cc{kk+1}+cc{kk+2})/3;
            mm=mm+1;
        end
    end
    for nn = 1:nrep
        kk = (k-1)*nrep +nn;
        if ismember('C',info{k,2}(:))
            TTV_c(j,:) = cc{kk};
            j = j+1;
        else
            if ismember('A',info(k,2))
                TTV_a(o,:) = cc{kk};
                o = o+1;
            elseif ismember('M',info(k,2)) 
                TTV_m(l,:) = cc{kk};
                l = l+1;
            elseif ismember('P',info(k,2))
                TTV_p(m,:) = cc{kk};
                m = m+1;
            end
        end
    end
 end
  %% Heatmap plot total volumne
 %tt = 135;
 tt=45;
bintime = 1;
stim = 20;
nbin = stim/bintime;
st =(-nbin+0.5)*bintime;
et =tt-(nbin+0.5)*bintime;
% st =(-nbin+0.5)*bintime;
% et =tt-(nbin+0.5)*bintime;
time = st:bintime:et;
figure;
clims =[-10,10]

colormap('jet')
tcl=tiledlayout(3,1);
nexttile
imagesc(TTV_a(1:end,40:85),'XData',time,clims);
title('anterior');
nexttile
%imagesc(TTV_m,'XData',time,clims);
%title('intermediate');
%nexttile
imagesc(TTV_p(1:end,40:85),'XData',time,clims);
title('posterior');
% nexttile
% imagesc(TTV_c(1:end,40:85),'XData',time,clims);
% title('control');

cb = colorbar;
cb.Layout.Tile = 'east';
title(tcl,' MAP deltaBP single trial, bin=5s, stimulation 0-60s ', fontsize=12);
%% Heatmap plot total volumne_averaged
 % tt = 135;
tt=35;
bintime = 1;
stim = 25; %stim=60;
nbin = stim/bintime
%st =(-nbin+0.5)*bintime;
st =(-10+0.5)*bintime;
et =tt-(10+0.5)*bintime;
time = st:bintime:et
figure;
clims =[-10,10];
showperiod = [51:85]

colormap('jet')
tcl=tiledlayout(4,1);
nexttile
imagesc(mTTV_a([1:5 7:8],showperiod ),'XData',time,clims);
title('anterior');
nexttile
imagesc(mTTV_m([1:5 7:10],showperiod ),'XData',time,clims);
title('intermediate');
nexttile
imagesc(mTTV_p([1:7 9:13 15:19],showperiod ),'XData',time,clims);
title('posterior');
nexttile
imagesc(mTTV_c(1:end,showperiod),'XData',time,clims);
title('control');

cb = colorbar;
cb.Layout.Tile = 'east';
title(tcl,' delta BP mouse average 3 trials, bin=5s, stimulation 0-25s ', fontsize=12);


%% Heatmap plot total volumne_averaged _ RANKORDERED/Sort
 % tt = 135;
tt=40;
bintime = 1;
stim = 30; %stim=60;
nbin = stim/bintime
%st =(-nbin+0.5)*bintime;
st =(-10+0.5)*bintime;
et =tt-(10+0.5)*bintime;
time = st:bintime:et
figure;
clims =[-10,10];
showperiod = [51:90]

dataA=mTTV_a([1:5 7:8],:)
% Rank rows by their mean values (descending)
rowMeans = mean(dataA(:,61:85), 2);
[~, rowOrder] = sort(rowMeans, 'descend');
data_sorted_A = dataA(rowOrder, :);


dataM=mTTV_m([1:5 7:10],:)
% Rank rows by their mean values (descending)
rowMeans = mean(dataM(:,61:85), 2);
[~, rowOrder] = sort(rowMeans, 'descend');
data_sorted_M = dataM(rowOrder, :);

dataP=mTTV_p([1:7 9:13 15:19],:)
% Rank rows by their mean values (descending)
rowMeans = mean(dataP(:,61:85), 2);
[~, rowOrder] = sort(rowMeans, 'descend');
data_sorted_P = dataP(rowOrder, :);


dataC=mTTV_c(1:end,:)
% Rank rows by their mean values (descending)
rowMeans = mean(dataC(:,61:85), 2);
[~, rowOrder] = sort(rowMeans, 'descend');
data_sorted_C = dataC(rowOrder, :);

colormap('jet')
tcl=tiledlayout(4,1);
nexttile
imagesc(data_sorted_A(:,showperiod ),'XData',time,clims);
title('anterior');
nexttile
imagesc(data_sorted_M(:,showperiod ),'XData',time,clims);
title('intermediate');
nexttile
imagesc(data_sorted_P(:,showperiod ),'XData',time,clims);
title('posterior');
nexttile
imagesc(data_sorted_C(:,showperiod),'XData',time,clims);
title('control');

cb = colorbar;
cb.Layout.Tile = 'east';
title(tcl,' all 4 group delta BP mouse average 3 trials, bin=5s, stimulation 0-25s ', fontsize=12);
%% Heatmap plot total volumne_averaged _ RANKORDERED/Sort
 % tt = 135;
tt=45;
bintime = 1;
stim = 25; %stim=60;
nbin = stim/bintime
%st =(-nbin+0.5)*bintime;
st =(-20+0.5)*bintime;
et =tt-(20+0.5)*bintime;
time = st:bintime:et
figure;
clims =[-10,10];
showperiod = [41:85]

dataA=mTTV_a([1:5 7:8],:)
% Rank rows by their mean values (descending)
rowMeans = mean(dataA(:,61:85), 2);
[~, rowOrder] = sort(rowMeans, 'descend');
data_sorted_A = dataA(rowOrder, :);


dataP=mTTV_p([1:7 9:13 15:19],:)
% Rank rows by their mean values (descending)
rowMeans = mean(dataP(:,61:85), 2);
[~, rowOrder] = sort(rowMeans, 'descend');
data_sorted_P = dataP(rowOrder, :);


colormap('jet')
tcl=tiledlayout(2,1);
nexttile
imagesc(data_sorted_A(:,showperiod ),'XData',time,clims);
title('anterior');
nexttile
% imagesc(mTTV_m([1:5 7:10],showperiod ),'XData',time,clims);
% title('intermediate');
% nexttile
imagesc(data_sorted_P(:,showperiod ),'XData',time,clims);
title('posterior');
% nexttile
% imagesc(mTTV_c(1:end,showperiod),'XData',time,clims);
% title('control');

cb = colorbar;
cb.Layout.Tile = 'east';
title(tcl,' delta BP mouse average 3 trials, bin=5s, stimulation 0-25s ', fontsize=12);

%% Plot Blood preesure (bin =1s=)
    tt = 40;
    bintime = 1;
    stim = 30;
    nbin = stim/bintime;
    st =(-10)*bintime;
    et =tt-(10)*bintime;
    time = st:bintime:et;
    showperiod = [50:90]
    % ylim([lowbound,topbound])
    % xlim([(-nbin+1)*bintime (tt/bin-nbin)*bintime])
    % set(gca,'FontSize', 15);
    % hold on
    % 
    % fill([0,0,120,120],[lowbound,topbound,topbound,lowbound],'yellow','FaceAlpha',0.2,'linestyle','none');
    % 
    lowbound = -15
    topbound = 7
    exchange_list='0123456789ABCDEF#'
    %colors ={'#0072BD','#D95319','#EDB120',	'#7E2F8E',	'#77AC30','#4DBEEE','#A2142F','#FF69B4'}
    colors = { '#E274A9', '#13AF68','#0F80FF','#424242'};  % line colors %"#E274A9" "#26A7E1'#FFE009'

    x=time
    h = figure('visible','on');
    hA = axes(h);

    mm=[]
    error=[]
    var = data_sorted_A(:,showperiod );
    mm(:,1) = mean(var,'omitnan')
    error(:,1)= std(var,'omitnan')/sqrt(size(var,1));
    var = data_sorted_M(:,showperiod );
    mm(:,2) = mean(var,'omitnan')
    error(:,2)= std(var,'omitnan')/sqrt(size(var,1));
    var = data_sorted_P(:,showperiod );
    mm(:,3) = mean(var,'omitnan')
    error(:,3)= std(var,'omitnan')/sqrt(size(var,1));
    var = data_sorted_C(:,showperiod );
    mm(:,4) = mean(var,'omitnan')
    error(:,4)= std(var,'omitnan')/sqrt(size(var,1));
    count=1
    p=[]
    for i=1:size(mm,2)
        y1=mm(:,i)+error(:,i)
        y2=mm(:,i)-error(:,i)
        x2=[x,fliplr(x)];
        inBetween = [y1',fliplr(y2')];
        string=(colors{count})
        num=zeros(1,3)
        for r=1:3
            temp_Coe1=find(exchange_list==string(r*2))-1;
            temp_Coe2=find(exchange_list==string(r*2+1))-1;
            num(r)=(16*temp_Coe1+temp_Coe2)/255;
        end
        fi=fill(x2,inBetween,num,'HandleVisibility','off');
        set(fi,'edgealpha',0,'facealpha',0.3)
        hold on
        plot(x,mm(:,i), 'LineWidth',1, 'Color', colors{count});

        count=count+1
        ylim([lowbound,topbound])
        xlim([st et])
        set(gca,'FontSize', 15);
    end
    hold on 
    ylims = ylim;  % Get current y-axis limits
    line([0 0], ylims, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1.5);
    xlims = xlim;  % Get current y-axis limits
    hold on
    line(xlims,[0 0],  'Color', 'k', 'LineStyle', '--', 'LineWidth', 1.5);
   
    group = {'Anterior','intermediate','posterior','Control'}
    lgd=legend(group,'Location','NorthEastOutside','Box','off')
    title(lgd,'Group')
    
    %fill([0,0,stim,stim],[lowbound,topbound,topbound,lowbound],'yellow','FaceAlpha',0.2,'linestyle','none');


%% Plot gut pressure raw trace 120s
time = -59.9998:0.0002:150 
h = figure('visible','on');
hA = axes(h);

for i = 1:15
    plot(time,data_extract{7,3}(:,i)+i*50,'color',"#A2142F",'LineWidth',1)
    hold on
end
for i = 16:30
    plot(time,data_extract{7,3}(:,i)+i*50,'color','#77AC30','LineWidth',1)
    hold on
end
for i = 31:45
    plot(time,data_extract{7,3}(:,i)+i*50,'color',"#0072BD",'LineWidth',1)
    hold on
end

lowbound = -0.5*50
topbound = 45*50
fill([0,0,120,120],[lowbound,topbound,topbound,lowbound],'blue','FaceAlpha',0.05,'linestyle','none');
hold on
    %fill([10,10,11,11],[lowbound,topbound,topbound,lowbound],'yellow','FaceAlpha',0.2,'linestyle','none');
    %hold on 
j = 0

xlim([-59.9998,150])
ylim([lowbound,topbound])
xlabel('Time(s)');
set(gca,'FontSize', 15);
    
    %set(hA, 'XTick', [], 'XTickLabel', []);
set(gcf,'units','points','position',[50,50,600,700])
set(hA, 'YTick', [], 'YTickLabel', []);
set(get(hA, 'XAxis'), 'Visible', 'on');
set(get(hA, 'YAxis'), 'Visible', 'off');
saveas(h,strcat( '120s_5Hz_gutpressure.pdf')); 

%% Plot Breathing raw trace 10s
time = -4.9998:0.0002:12.5 
h = figure('visible','on');
hA = axes(h);

group = {'anterior','intermediate','posterior'}
lgd=legend(group,'Location','NorthEastOutside','Box','off')
title(lgd,'Group')


for i = 1:15
    plot(time,data_extract{3,4}(:,i)+i*0.5,'color',"#A2142F",'LineWidth',1)
    hold on
end
for i = 16:30
    plot(time,data_extract{3,4}(:,i)+i*0.5,'color','#77AC30','LineWidth',1)
    hold on
end
for i = 31:45
    plot(time,data_extract{3,4}(:,i)+i*0.5,'color',"#0072BD",'LineWidth',1)
    hold on
end

lowbound = -0.5
topbound = 23
fill([0,0,10,10],[lowbound,topbound,topbound,lowbound],'blue','FaceAlpha',0.05,'linestyle','none');
hold on
    %fill([10,10,11,11],[lowbound,topbound,topbound,lowbound],'yellow','FaceAlpha',0.2,'linestyle','none');
    %hold on 
j = 0

xlim([-4.9998,12.5])
ylim([lowbound,topbound])
xlabel('Time(s)');
set(gca,'FontSize', 15);
    
    %set(hA, 'XTick', [], 'XTickLabel', []);
set(gcf,'units','points','position',[50,50,600,700])
set(hA, 'YTick', [], 'YTickLabel', []);
set(get(hA, 'XAxis'), 'Visible', 'on');
set(get(hA, 'YAxis'), 'Visible', 'off');
saveas(h,strcat( 'BreathTrace.pdf')); 



%% Plot Breathing raw trace 120s 5Hz
time = -59.9998:0.0002:150 
h = figure('visible','on');
hA = axes(h);


for i = 1:15
    plot(time,data_extract{7,4}(:,i)+i,'color',"#A2142F",'LineWidth',1)
    hold on
end
for i = 16:30
    plot(time,data_extract{7,4}(:,i)+i,'color','#77AC30','LineWidth',1)
    hold on
end
for i = 31:45
    plot(time,data_extract{7,4}(:,i)+i,'color',"#0072BD",'LineWidth',1)
    hold on
end

lowbound = -0.5
topbound = 45
fill([0,0,120,120],[lowbound,topbound,topbound,lowbound],'blue','FaceAlpha',0.05,'linestyle','none');
hold on
    %fill([10,10,11,11],[lowbound,topbound,topbound,lowbound],'yellow','FaceAlpha',0.2,'linestyle','none');
    %hold on 
j = 0

xlim([-59.9998,150])
ylim([lowbound,topbound])
xlabel('Time(s)');
set(gca,'FontSize', 15);
    
    %set(hA, 'XTick', [], 'XTickLabel', []);
set(gcf,'units','points','position',[50,50,600,700])
set(hA, 'YTick', [], 'YTickLabel', []);
set(get(hA, 'XAxis'), 'Visible', 'on');
set(get(hA, 'YAxis'), 'Visible', 'off');
saveas(h,strcat( '120s_5Hz_BreathTrace.pdf')); 

%% Plot Breathing raw trace 60s 20Hz
time = -29.9998:0.0002:75 
h = figure('visible','on');
hA = axes(h);


for i = 1:15
    plot(time,data_extract{2,4}(:,i)+i,'color',"#A2142F",'LineWidth',1)
    hold on
end
for i = 16:30
    plot(time,data_extract{2,4}(:,i)+i,'color','#77AC30','LineWidth',1)
    hold on
end
for i = 31:45
    plot(time,data_extract{2,4}(:,i)+i,'color',"#0072BD",'LineWidth',1)
    hold on
end

lowbound = -0.5
topbound = 45
fill([0,0,60,60],[lowbound,topbound,topbound,lowbound],'blue','FaceAlpha',0.05,'linestyle','none');
hold on
    %fill([10,10,11,11],[lowbound,topbound,topbound,lowbound],'yellow','FaceAlpha',0.2,'linestyle','none');
    %hold on 
j = 0

xlim([-29.9998,75])
ylim([lowbound,topbound])
xlabel('Time(s)');
set(gca,'FontSize', 15);
    
    %set(hA, 'XTick', [], 'XTickLabel', []);
set(gcf,'units','points','position',[50,50,600,700])
set(hA, 'YTick', [], 'YTickLabel', []);
set(get(hA, 'XAxis'), 'Visible', 'on');
set(get(hA, 'YAxis'), 'Visible', 'off');
saveas(h,strcat( '60s_20Hz_BreathTrace.pdf')); 

%% Plot Blood pressure raw trace 60s 20Hz
time = -59.9998:0.0002:75 
h = figure('visible','on');
hA = axes(h);


for i = 1:15
    plot(time,data_extract{8,2}(:,i)+i*150,'color',"#A2142F",'LineWidth',1)
    hold on
end

lowbound = -0.5
topbound = 150*15
fill([0,0,60,60],[lowbound,topbound,topbound,lowbound],'blue','FaceAlpha',0.05,'linestyle','none');
hold on
    %fill([10,10,11,11],[lowbound,topbound,topbound,lowbound],'yellow','FaceAlpha',0.2,'linestyle','none');
    %hold on 
j = 0

xlim([-59.9998,75])
ylim([lowbound,topbound])
xlabel('Time(s)');
set(gca,'FontSize', 15);
    
    %set(hA, 'XTick', [], 'XTickLabel', []);
set(gcf,'units','points','position',[50,50,600,700])
set(hA, 'YTick', [], 'YTickLabel', []);
set(get(hA, 'XAxis'), 'Visible', 'on');
set(get(hA, 'YAxis'), 'Visible', 'off');
saveas(h,strcat( '60s_20Hz_BP.pdf')); 

figure
for i = 16:30
    plot(time,data_extract{8,2}(:,i)+i*150,'color','#77AC30','LineWidth',1)
    hold on
end
lowbound = -0.5
topbound = 150*15
fill([0,0,60,60],[lowbound,topbound,topbound,lowbound],'blue','FaceAlpha',0.05,'linestyle','none');
hold on
    %fill([10,10,11,11],[lowbound,topbound,topbound,lowbound],'yellow','FaceAlpha',0.2,'linestyle','none');
    %hold on 
j = 0

xlim([-59.9998,75])
ylim([lowbound,topbound])
xlabel('Time(s)');
set(gca,'FontSize', 15);
    
    %set(hA, 'XTick', [], 'XTickLabel', []);
set(gcf,'units','points','position',[50,50,600,700])
set(hA, 'YTick', [], 'YTickLabel', []);
set(get(hA, 'XAxis'), 'Visible', 'on');
set(get(hA, 'YAxis'), 'Visible', 'off');
saveas(h,strcat( '60s_20Hz_BP.pdf')); 

figure
for i = 31:45
    plot(time,data_extract{8,2}(:,i)+i*150,'color',"#0072BD",'LineWidth',1)
    hold on
end

lowbound = -0.5
topbound = 150*15
fill([0,0,60,60],[lowbound,topbound,topbound,lowbound],'blue','FaceAlpha',0.05,'linestyle','none');
hold on
    %fill([10,10,11,11],[lowbound,topbound,topbound,lowbound],'yellow','FaceAlpha',0.2,'linestyle','none');
    %hold on 
j = 0

xlim([-59.9998,75])
ylim([lowbound,topbound])
xlabel('Time(s)');
set(gca,'FontSize', 15);

%set(hA, 'XTick', [], 'XTickLabel', []);

set(gcf,'units','points','position',[50,50,600,700])
set(hA, 'YTick', [], 'YTickLabel', []);
set(get(hA, 'XAxis'), 'Visible', 'on');
set(get(hA, 'YAxis'), 'Visible', 'off');
saveas(h,strcat( '60s_20Hz_BP.pdf')); 
%%

% for i = 1:15   
%     figure
%     plot(SBP_zscore{i})
%     %plot(bandstop(data_extract{2,2}(:,i),[0.2 3],sampleRate));
% end
%88:132 batch3
for i = 13:15
    %mm = (i-1)*3+1
    figure
    %plot(data_extract{1,4}(:,i))
    plot(detrend(data_extract{1,4}(:,i)))
end


%% Final plot raw trace_BP

j=0
o=0
l=0
m=0

figure
tcl=tiledlayout(4,15);
c=3
for i = 1:44
    mm = (i-1)*3+c
    if ismember('C',info{i,2}(:))
        j = j+1;
        nexttile(3*15+j);
    else
        if ismember('A',info(i,2))
            o = o+1;
            nexttile(o);
        elseif ismember('M',info(i,2)) 
            l = l+1;
            nexttile(15+l);
        elseif ismember('P',info(i,2))
            m = m+1;
            nexttile(30+m);
        end
    end
    plot(data_extract{2,1}(:,mm));
    %plot(detrend(data_extract{8,4}(:,mm)));
     %plot(detrend(data_extract{2,4}(250000:350000,mm)));
     title([info{i,2},' M',info{i,3},' T',num2str(c)]);
     %plot(data_extract{1,4}(:,i))

end
title(tcl,' BP-60s', fontsize=12);
%% Final plot raw trace_BR

j=0
o=0
l=0
m=0

figure
tcl=tiledlayout(4,15);
c=3
for i = 1:44
    mm = (i-1)*3+c
    if ismember('C',info{i,2}(:))
        j = j+1;
        nexttile(3*15+j);
    else
        if ismember('A',info(i,2))
            o = o+1;
            nexttile(o);
        elseif ismember('M',info(i,2)) 
            l = l+1;
            nexttile(15+l);
        elseif ismember('P',info(i,2))
            m = m+1;
            nexttile(30+m);
        end
    end
    plot(data_extract{1,4}(:,mm));
    %plot(detrend(data_extract{8,4}(:,mm)));
     %plot(detrend(data_extract{2,4}(250000:350000,mm)));
     title([info{i,2},' M',info{i,3},' T',num2str(c)]);
     %plot(data_extract{1,4}(:,i))

end
title(tcl,' BR-10s', fontsize=12);
%% BP- Final plot raw trace_BP_in one graph

j=0
o=0
l=0
m=0
%time = -59.9998:0.0002:75
time = -9.9998:0.0002:70; 
snn= 50.0002*sampleRate;
enn = 130*sampleRate;
showwindow = int32([snn:1:enn]);
figure
ii=0
tracelabel={}
for i = 1:44
    mm = (i-1)*3;
    if ismember('P',info(i,2))
    %if ismember('C',info{i,2}(:))
        ii=ii+1
        %if (ii ~=6)   
         if (ii ~=8) &(ii~=14)
            for jj =1:3
                  plot(time,data_extract{2,2}(showwindow,mm+jj)+o*30);
                  o = o+1;
                  tracelabel{o} = [num2str(info{i,1}),'-',info{i,2},'-',info{i,3},'-T',num2str(jj)]
                  hold on
            end
        end
    end
end
legend(tracelabel);
xregion(0, 60);

plot([-10 0],[50 50],'k-','linewidth',3); % 3point for the bar
text(-5,50,'10s','horiz','center','vert','top'); 

% plot([-1 -1 0],[2.1 2 2],'k-','linewidth',3); % 3point for the bar
% text(-0.5,1.5,'1s','horiz','center','vert','top'); 
% text(80,.85,'0.1 y units ','horiz','right','vert','middle');
fname ='BP-1min-raw'
title(fname, fontsize=12);

%% Breathing- Final plot raw trace_BR_in one graph

j=0
o=0
l=0
m=0
time = -9.9998:0.0002:12.5 
figure
ii=0
tracelabel={}
for i = 1:44
    mm = (i-1)*3;
    %if ismember('P',info(i,2))
    if ismember('C',info{i,2}(:))
        ii=ii+1
        %if (ii ~=6)   
        % if (ii ~=8) &(ii~=14)
            for jj =1:3
                  plot(time,data_extract{1,4}(:,mm+jj)-o*1);
                  o = o+1;
                  tracelabel{o} = [num2str(info{i,1}),'-',info{i,2},'-',info{i,3},'-T',num2str(jj)]
                  hold on
            end
        %end
    end
end
legend(tracelabel);
xregion(0, 10);

plot([-1 0],[2 2],'k-','linewidth',3); % 3point for the bar
text(-0.5,1.5,'1s','horiz','center','vert','top'); 

% plot([-1 -1 0],[2.1 2 2],'k-','linewidth',3); % 3point for the bar
% text(-0.5,1.5,'1s','horiz','center','vert','top'); 
% text(80,.85,'0.1 y units ','horiz','right','vert','middle');
fname ='BR-10s-raw'
title(fname, fontsize=12);


%% Breathing-Select 8 representative trace
[fileDir] = uigetdir( 'Output');
 cd(fileDir);
j=0
o=0
l=0
m=0
time = -4.9998:0.0002:12.5; 
snn= 5.0002*sampleRate;
enn = 22.5*sampleRate;
showwindow = int32([snn:1:enn]);
figure
colors = { '#E274A9','#E274A9','#E274A9','#E274A9','#E274A9','#E274A9', '#13AF68','#13AF68', '#13AF68','#13AF68', '#13AF68','#13AF68','#0F80FF','#0F80FF','#0F80FF','#0F80FF','#0F80FF','#0F80FF','#424242','#424242','#424242','#424242','#424242','#424242'};  % line colors %"#E274A9" "#26A7E1'#FFE009'

ii=0
tracelabel={}
%ii = [2 18 22 29]
%jj = [1 2 2 1]
ii = [2 2 4 4 5 5 18 18 10 10 6 6 45 45 22 22 11 11 25 25 28 28 26 26]
jj = [1 2 2 3 2 3 2 3 2 3 2 3 2 3 1 2 2 3 2 3 2 3 1 2]

for i = 1:24
    mm = (ii(i)-1)*3;
    plot(time,data_extract{1,4}(showwindow,mm+jj(i))-o*1,'Color',colors{i});
    o = o+1;
    tracelabel{o} = [info{ii(i),2},'-',info{ii(i),3},'-T',num2str(jj(i))]
    hold on           
end

legend(tracelabel);
xregion(0, 10);
ylim([-25,3])
plot([-1 0],[2 2],'k-','linewidth',3); % 3point for the bar
text(-0.5,1.5,'1s','horiz','center','vert','top'); 

% plot([-1 -1 0],[2.1 2 2],'k-','linewidth',3); % 3point for the bar
% text(-0.5,1.5,'1s','horiz','center','vert','top'); 
% text(80,.85,'0.1 y units ','horiz','right','vert','middle');
fname ='BR-10s-example trace'
title(fname, fontsize=12);

print(gcf,'-vector','-dsvg',[fname,'.svg']) % svg
%print(gcf,'-vector','-dpdf',[fname,'.pdf']) % pdf

 %%
fs = 1/SampleRate*1000

BR_base=[]
BR = data_extract{1,4}(:,1)
BR_base = bandpass(BR,[0.1 0.3],fs)
plot(BR)
hold on
plot(BR_base)

%% BP-Select 2*2*3 representative trace
[fileDir] = uigetdir( 'Output');
 cd(fileDir);
j=0
o=0
l=0
m=0
downsample=500 %downsample times
time = -29.9998:0.0002* downsample:60; 
snn= (30+0.0002)*(sampleRate);
enn = 120*(sampleRate);
% time = -14.9998:0.0002* downsample:30; 
% snn= (45+0.0002)*(sampleRate);
% enn = 90*(sampleRate);
showwindow = int32([snn:1:enn]);
figure
%colors = { '#E274A9','#E274A9', '#E274A9','#E274A9','#E274A9','#E274A9','#13AF68','#13AF68','#0F80FF','#0F80FF','#424242','#424242'};  % line colors %"#E274A9" "#26A7E1'#FFE009'
%colors =  {"#424242","#424242","#424242","#424242","#0000FF","#0000FF","#0000FF","#0000FF","#FF0000","#FF0000","#FF0000","#FF0000","#FF0000"};
colors = { '#E274A9','#E274A9','#E274A9','#E274A9','#E274A9','#E274A9', '#13AF68','#13AF68', '#13AF68','#13AF68', '#13AF68','#13AF68','#0F80FF','#0F80FF','#0F80FF','#0F80FF','#0F80FF','#0F80FF','#424242','#424242','#424242','#424242','#424242','#424242'};  % line colors %"#E274A9" "#26A7E1'#FFE009'

ii=0
tracelabel={}
%ii = [2 18 22 29]
%jj = [1 2 2 1]
ii = [2 2 3 3 4 4 6 6 10 10 34 34 19 19 37 37 38 38 43 43 27 27 25 25]
jj = [2 3 2 3 2 3 2 3 2 3 2 3 1 3 2 3 2 3 1 2 2 3 1 2 ]

cd = data_extract{2,2};

for i = 1:length(ii)
    mm = (ii(i)-1)*3;
    current=cd(showwindow,mm+jj(i));
    y0 = smoothdata(current,'movmean',50);
    y = decimate(y0, downsample);
    plot(time,y-l-max(y),'Color',colors{i});
    l = l+(max(y)-min(y))*1.1;
    o = o+1;
    tracelabel{o} = [info{ii(i),2},'-',info{ii(i),3},'-T',num2str(jj(i))];
    hold on
           
end
legend(tracelabel);
xregion(0, 60);

plot([-5 -5], [2 4],'k-',[-5 0],[2 2],'k-','linewidth',3); % 3point for the bar
text(-4.5,1.5,'5s','horiz','center','vert','top'); 
text(-5.5,4.5,'2cmH2O','horiz','right','vert','top')
% plot([-1 -1 0],[2.1 2 2],'k-','linewidth',3); % 3point for the bar
% text(-0.5,1.5,'1s','horiz','center','vert','top'); 
% text(80,.85,'0.1 y units ','horiz','right','vert','middle');
fname ='BP-1min-all-example trace'
title(fname, fontsize=12);

print(gcf,'-vector','-dsvg',[fname,'.svg']) % svg
%print(gcf,'-vector','-dpdf',[fname,'.pdf']) % pdf


%% plot BP zscore
nrep=3;
j=1;
o=1;
l=1;
m=1;
jj=1;
oo=1;
ll=1;
mm=1;
TTV_a = [];
TTV_m = [];
TTV_p = [];
TTV_c = [];
mTTV_a = [];
mTTV_m = [];
mTTV_p = [];
mTTV_c = [];

cc = MAP_zscore;
for k = 1:size(info,1)
    kk = (k-1)*nrep +1;
    cc_zscore(k,1)=mean(cc{kk}(76:85));
    cc_zscore(k,2)=mean(cc{kk+1}(76:85));
    cc_zscore(k,3)=mean(cc{kk+2}(76:85));
end
% 
%     if ismember('C',info{k,2}(:))
%         mTTV_c(jj,:) = (cc{kk}+cc{kk+1}+cc{kk+2})/3;
%         jj=jj+1;
%     else
%         if ismember('A',info(k,2))
%             mTTV_a(oo,:) =   (cc{kk}+cc{kk+1}+cc{kk+2})/3;
%             oo=oo+1;
%         elseif ismember('M',info(k,2)) 
%             mTTV_m(ll,:) =  (cc{kk}+cc{kk+1}+cc{kk+2})/3;
%             ll=ll+1;
%         elseif ismember('P',info(k,2))
%             mTTV_p(mm,:) =   (cc{kk}+cc{kk+1}+cc{kk+2})/3;
%             mm=mm+1;
%         end
%     end
%     for nn = 1:nrep
%         kk = (k-1)*nrep +nn;
%         if ismember('C',info{k,2}(:))
%             TTV_c(j,:) = cc{kk};
%             j = j+1;
%         else
%             if ismember('A',info(k,2))
%                 TTV_a(o,:) = cc{kk};
%                 o = o+1;
%             elseif ismember('M',info(k,2)) 
%                 TTV_m(l,:) = cc{kk};
%                 l = l+1;
%             elseif ismember('P',info(k,2))
%                 TTV_p(m,:) = cc{kk};
%                 m = m+1;
%             end
%         end
%     end
% end 



