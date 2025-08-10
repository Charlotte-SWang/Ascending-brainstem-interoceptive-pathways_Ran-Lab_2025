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
%cc=BR_CI;
cc =CI_amp_tv;
%cc=amp_stage_raw;
%% plot TV zscore
nrep=3
j=1;
o=1;
l=1;
m=1;
jj=1;
oo=1;
ll=1;
mm=1;

TTV_v = [];
TTV_p = [];
TTV_c = [];

mTTV_c = [];
mTTV_p = [];
mTTV_v = [];
cc=[]
cc = single_ci_amp_tv


info_current = trialinfo_extract{3,2};
[unique_vals, ~, idx] = unique([info_current{:,1}],'stable');
counts = accumarray(idx, 1);  % returns 5


for k = 1:size(counts,1)
    ccsubset =cc((sum(counts(1:k))-counts(k)+1):sum(counts(1:k)))
    ccmat= cat(1,ccsubset{:})

    for nn = 1:counts(k)
        if k==1
            kk = nn;
        else
            kk = sum(counts(1:k-1)) +nn;
        end
        if ismember('C',info_current{kk,4})
            TTV_c(j,:) = cc{kk};
            j = j+1;
        else
            if ismember('PVH',info_current{kk,3})
                TTV_p(o,:) = cc{kk};
                o = o+1;
            elseif ismember('VLM',info_current{kk,3})
                TTV_v(l,:) = cc{kk};
                l = l+1;
            end
        end
    end
    if ismember('C',info_current{kk,4})
        mTTV_c(jj,:) = sum(ccmat,1,"omitnan")/counts(k);
        jj=jj+1;
    else
        if ismember('PVH',info_current{kk,3})
            mTTV_p(oo,:) =  sum(ccmat,1,"omitnan")/counts(k);
            oo=oo+1;
        elseif ismember('VLM',info_current{kk,3})
            mTTV_v(ll,:) = sum(ccmat,1,"omitnan")/counts(k);
            ll=ll+1;
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
clims =[-0.2,0.2];
%clims = [-1,1]
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
TTV_v = [];
TTV_p = [];
TTV_c = [];

mTTV_c = [];
mTTV_p = [];
mTTV_v = [];
gmTTV={};

gmTTV{1}=[];
gmTTV{2}=[];
gmTTV{3}=[];

groupinfo=[];
imTTV={};

cc =single_ci_amp_tv;

info_current = trialinfo_extract{3,2};
[unique_vals, ~, idx] = unique([info_current{:,1}],'stable');
counts = accumarray(idx, 1);  % returns 5


for k = 1:size(counts,1)
    ccmat = {};
  
        if k==1
            kk = 1;
        else
            kk = sum(counts(1:k-1)) +1;
        end
        ccmat = [cc{kk}; cc{kk+1}; cc{kk+2}];
 
    if ismember('C',info_current{sum(counts(1:k)),4})
        mTTV_c{jj} =ccmat;
        gmTTV{1}=[gmTTV{1}; mTTV_c{jj}];
       
        imTTV=[imTTV; mTTV_c{jj}];
        groupinfo=[groupinfo 1];

        jj=jj+1;
        j=j+1;

    else
        if ismember('PVH',info_current{sum(counts(1:k)),3})
            mTTV_p{oo} = ccmat;
            gmTTV{2}=[gmTTV{2};mTTV_p{oo}];

            imTTV=[imTTV; mTTV_p{oo}];

            groupinfo=[groupinfo 2];

            oo=oo+1;
            o = o+1;
        elseif ismember('VLM',info_current{sum(counts(1:k)),3})
            mTTV_v{mm} = ccmat;
            gmTTV{3}=[gmTTV{3}; mTTV_v{mm}];

            imTTV=[imTTV; mTTV_v{mm}];
            groupinfo=[groupinfo 3];

            mm=mm+1;
        end   
    end
end


%% plot BREATHING tidal volume-group 3d histogram
figure
xs = 0:5:10;
colors = jet(numel(xs));
for i = 1:numel(xs)
    histogram2(repmat(xs(i),1,size(gmTTV{i},1)),gmTTV{i}.','FaceColor',colors(i,:),'Normalization','percentage')
    hold on
end
hold off

%% plot BREATHING tidal volume- individual 3d histogram
figure
xs = 1:1:size(imTTV,1)
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
xs = 1:1:size(imTTV,1)
colormap('jet')
imTTV_group={}
[groupinfo_sorted, rowOrder] = sort(groupinfo, 'ascend');
imTTV_sorted = imTTV(rowOrder);

for i = 1:3
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
subplot(3, 1, 1)
histogram(gmTTV{1}.','BinWidth',0.01,'Normalization','count','FaceColor','Black')
xlim([-0.2 0.5])
ylim([0 80])
subplot(3, 1, 2)
histogram(gmTTV{2}.','BinWidth',0.01,'Normalization','count','FaceColor','Blue')
xlim([-0.2 0.5])
ylim([0 150])
subplot(3, 1, 3)
histogram(gmTTV{3}.','BinWidth',0.01,'Normalization','count','FaceColor','Red')
xlim([-0.2 0.5])
ylim([0 30])
%caxis([0 100])
hold off

%% plot BREATHING tidal volume- Breathing Tidal Volume-Preview of Kernel Fitting
figure 
% Plot normalized histogram (as a density, not counts)
normaldata = gmTTV{1};
normalpd = fitdist(normaldata,"Kernel");
plot(normalpd);
%% plot BREATHING tidal volume- Breathing Tidal Volume-Final Cumulative distribution
%colors = { "#E274A9", '#13AF68',"#26A7E1","#E95412"};  % line colors %"#E274A9" "#26A7E1'#FFE009'
%colors = { "#E274A9", '#13AF68',"#0F80FF","#424242"};  % line colors %"#E274A9" "#26A7E1'#FFE009'
colors = {"#424242","#0000FF","#FF0000"};
hold on;

% Define common x-range
%x = linspace(min(cellfun(@min, gmTTV)), max(cellfun(@max, gmTTV)), 200);
x = linspace(-1,1,200)
for i = 1:3
    data = gmTTV{i};
    pd = fitdist(data, 'Kernel');
    y = cdf(pd, x);  % evaluate CDF
    plot(x, y, 'Color', colors{i}, 'LineWidth',2);  % clean line only
end
xlim([-0.7,0.7])
xlabel('Change index');
ylabel('Cumulative distribution');
legend('Control', 'PVH', 'VLM');
title('Overlaid CDFs from Kernel Distributions of 5s tidal volume data');
grid off;
%% plot BREATHING tidal volume-Breathing Tidal Volume-Density Cumulative distribution
figure
% Plot kernel density estimate
%colors = { "#E274A9", '#13AF68',"#0F80FF","#424242"};  % line colors %"#E274A9" "#26A7E1'#FFE009'
%colors = { "#E95412", '#13AF68',"#26A7E1",'#FFE009'};  % line colors %"#E274A9" "#26A7E1
colors = {"#424242","#0000FF","#FF0000"};
hold on;

% Define common x-range
%x = linspace(min(cellfun(@min, gmTTV)), max(cellfun(@max, gmTTV)), 200);
x = linspace(-1,0.6,200)
for i = 1:3
    [f, x] = ksdensity(gmTTV{i});
    plot(x, f, 'r-', 'LineWidth', 2,'Color',colors{i});
    hold on;
end

% Labels and formatting
xlabel('Change index');
ylabel('Density');
xlim([-0.6,0.6])
legend('Control', 'PVH', 'VLM');
title('Histogram with Kernel Density Estimate of 5s tidal volume data');

%% Normal KS test for 5s breathing tidal volume Pooled
g1=gmTTV{1}
g2=gmTTV{2}
g3=gmTTV{3}
[h2, p2, stat2] = kstest2(g1, g2, 'Tail', 'unequal')

[h3, p3, stat3] = kstest2(g1, g3, 'Tail', 'unequal')


%% bPermuted KS test for 5s breathing tidal volume--ADD FUNCTION
rng(1120)
[pathname_function] = uigetdir( 'Select the physiology function');
addpath(pathname_function)
%% Permuted KS test for 5s breathing tidal volume
data = [imTTV_group{1};imTTV_group{2}];
labels = [ones(length(imTTV_group{1}),1); 2*ones(length(imTTV_group{2}),1)];
[stat, p] = exact_ks_permutation_test(data,labels, 'two-sided',100000);
adjust_p = p*2

data = [imTTV_group{1};imTTV_group{3}];
labels = [ones(length(imTTV_group{1}),1); 2*ones(length(imTTV_group{3}),1)];
[stat, p] = exact_ks_permutation_test(data,labels, 'two-sided',100000);
adjust_p = p*2

%
data = [imTTV_group{1};imTTV_group{2}];
labels = [ones(length(imTTV_group{1}),1); 2*ones(length(imTTV_group{2}),1)];
[stat, p] = difference_permutation_test(data,labels, 'two-sided',100000);
adjust_p_diff = p*2

data = [imTTV_group{1};imTTV_group{3}];
labels = [ones(length(imTTV_group{1}),1); 2*ones(length(imTTV_group{3}),1)];
[stat, p] = difference_permutation_test(data,labels, 'two-sided',100000);
adjust_p_diff = p*2
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

%% Final plot raw trace_BR

j=0
o=0
l=0
m=0

figure
tcl=tiledlayout(3,16);
c=3
for i = 1:105
    if trialinfo_extract{3,2}{i,5} ~=c
        continue
    else
    %mm = (i-1)*3+c
        if ismember('C',trialinfo_extract{3,2}{i,4})
            j = j+1;
            nexttile(j);
        else
            if ismember('PVH',trialinfo_extract{3,2}{i,3})
                o = o+1;
                nexttile(16+o);
            elseif ismember('VLM',trialinfo_extract{3,2}{i,3})
                l = l+1;
                nexttile(2*16+l);
            end
        end
        plot(data_extract{3,2}(:,i));
        %plot(detrend(data_extract{8,4}(:,mm)));
         %plot(detrend(data_extract{2,4}(250000:350000,mm)));
         title([' M',trialinfo_extract{3,2}{i,2},trialinfo_extract{3,2}{i,3},trialinfo_extract{3,2}{i,4},' T',num2str(trialinfo_extract{3,2}{i,5})]);
         %plot(data_extract{1,4}(:,i))
    end
end
title(tcl,' BR-5s', fontsize=12);
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
time = -4.9998:0.001:10 
figure
ii=0
tracelabel={}
for i = 1:105
    mm = (i-1)*3;
    %if ismember('P',info(i,2))
    %if ismember('C',trialinfo_extract{3,2}{i,4})
    %if ismember('PVH',trialinfo_extract{3,2}{i,3}) & ismember('E',trialinfo_extract{3,2}{i,4})
    if ismember('VLM',trialinfo_extract{3,2}{i,3}) & ismember('E',trialinfo_extract{3,2}{i,4})
        %ii=ii+1
        %if (ii ~=6)   
        % if (ii ~=8) &(ii~=14)
           % for jj =1:3
                  plot(time,data_extract{3,2}(:,i)-o*0.5);
                  o = o+1;
                  tracelabel{o} = [num2str(trialinfo_extract{3,2}{i,1}),trialinfo_extract{3,2}{i,2},trialinfo_extract{3,2}{i,3},trialinfo_extract{3,2}{i,4},' T',num2str(trialinfo_extract{3,2}{i,5})]
                  hold on
            %end
        %end
    end
end
legend(tracelabel);
xregion(0, 5);

plot([-1 0],[2 2],'k-','linewidth',3); % 3point for the bar
text(-0.5,1.5,'1s','horiz','center','vert','top'); 

% plot([-1 -1 0],[2.1 2 2],'k-','linewidth',3); % 3point for the bar
% text(-0.5,1.5,'1s','horiz','center','vert','top'); 
% text(80,.85,'0.1 y units ','horiz','right','vert','middle');
fname ='BR-5s-raw-VLM'
title(fname, fontsize=12);

print(gcf,'-vector','-dsvg',[fname,'.svg']) % svg
%print(gcf,'-vector','-dpdf',[fname,'.pdf']) % pdf
%% Breathing-Select 8 representative trace
j=0
o=0
l=0
m=0
time = -4.999:0.001:10; 
snn= 0.001*sampleRate;
enn = 15*sampleRate;
showwindow = int32([snn:1:enn]);
figure
%colors = { '#E274A9','#E274A9', '#E274A9','#E274A9','#E274A9','#E274A9','#13AF68','#13AF68','#0F80FF','#0F80FF','#424242','#424242'};  % line colors %"#E274A9" "#26A7E1'#FFE009'
colors =  {"#424242","#424242","#424242","#424242","#424242","#424242","#0000FF","#0000FF","#0000FF","#0000FF","#0000FF","#0000FF","#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#FF0000"};
ii=0
tracelabel={}
%ii = [2 18 22 29]
%jj = [1 2 2 1]
ii = [52 54 86 87 2 3 103 104 20 21 22 24 37 38 73 74 32 33 ]
jj = [ 1 1 1 1 1 1 1 1 1 1 1 1 1 1 -1 -1 1 1]
cd = data_extract{3,2};
ci = trialinfo_extract{3,2};
for i = 1:length(ii)
    mm = ii(i);
    plot(time,cd(showwindow,mm)-o*1,'Color',colors{i});
    o = o+1;
    tracelabel{o} = [ci{mm,4},'-',ci{mm,3},'-',ci{mm,2},'-T',num2str(ci{mm,5}),'-idx',num2str(mm)];
    hold on
           

end
legend(tracelabel);
xregion(0, 5);

plot([-1 0],[2 2],'k-','linewidth',3); % 3point for the bar
text(-0.5,1.5,'1s','horiz','center','vert','top'); 

% plot([-1 -1 0],[2.1 2 2],'k-','linewidth',3); % 3point for the bar
% text(-0.5,1.5,'1s','horiz','center','vert','top'); 
% text(80,.85,'0.1 y units ','horiz','right','vert','middle');
fname ='BR-5s-example trace'
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


%% 5min gasp analysis -Select 2*2*3 representative trace
j=0
o=0
l=0
m=0
time = -299.999:0.001:315; 
snn= 0.001*sampleRate;
enn = 615*sampleRate;
showwindow = int32([snn:1:enn]);
figure
%colors = { '#E274A9','#E274A9', '#E274A9','#E274A9','#E274A9','#E274A9','#13AF68','#13AF68','#0F80FF','#0F80FF','#424242','#424242'};  % line colors %"#E274A9" "#26A7E1'#FFE009'
colors =  {"#424242","#424242","#424242","#424242","#0000FF","#0000FF","#0000FF","#0000FF","#FF0000","#FF0000","#FF0000","#FF0000","#FF0000"};
ii=0
tracelabel={}
%ii = [2 18 22 29]
%jj = [1 2 2 1]
ii = [2 3 84 85 102 103 64 65 37 38 71 72]

cd = data_extract{1,2};
ci = trialinfo_extract{1,2};
for i = 1:length(ii)
    mm = ii(i);
    plot(time,cd(showwindow,mm)-o*1,'Color',colors{i});
    o = o+1;
    tracelabel{o} = [ci{mm,4},'-',ci{mm,3},'-',ci{mm,2},'-T',num2str(ci{mm,5}),'-idx',num2str(mm)];
    hold on
           

end
legend(tracelabel);
xregion(0, 300);

plot([-5 -5], [2 2.19],'k-',[-30 0],[2 2],'k-','linewidth',3); % 3point for the bar
text(-4.5,1.5,'30s','horiz','center','vert','top'); 
text(-5.5,4.5,'NO UNIT','horiz','right','vert','top')
% plot([-1 -1 0],[2.1 2 2],'k-','linewidth',3); % 3point for the bar
% text(-0.5,1.5,'1s','horiz','center','vert','top'); 
% text(80,.85,'0.1 y units ','horiz','right','vert','middle');
fname ='BR-Gasp-5min-all-example trace'
title(fname, fontsize=12);

print(gcf,'-vector','-dsvg',[fname,'.svg']) % svg
%print(gcf,'-vector','-dpdf',[fname,'.pdf']) % pdf
