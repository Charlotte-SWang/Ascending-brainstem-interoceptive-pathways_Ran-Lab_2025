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
mTTV_v = [];
mTTV_m = [];
mTTV_p = [];
mTTV_c = [];
cc = Delta_BP


info_current = trialinfo_extract{2,3};
[unique_vals, ~, idx] = unique([info_current{:,1}],'stable');
counts = accumarray(idx, 1);  % returns 5


for k = 1:size(counts,1)
    ccsubset =cc((sum(counts(1:k))-counts(k)+1):sum(counts(1:k)))
    ccmat= cat(2,ccsubset{:})

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
        mTTV_c(jj,:) = sum(ccmat,2)/counts(k);
        jj=jj+1;
    else
        if ismember('PVH',info_current{kk,3})
            mTTV_p(oo,:) =  sum(ccmat,2)/counts(k);
            oo=oo+1;
        elseif ismember('VLM',info_current{kk,3})
            mTTV_v(ll,:) = sum(ccmat,2)/counts(k);
            ll=ll+1;
        end
    end
end


  %% Heatmap plot total volumne
 %tt = 135;
 tt=50;
 showperiod = 40:90

bintime = 1;
nstim = 20; % time before stimulation
nbin = nstim/bintime;
%st =(-nbin+0.5)*bintime;
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
imagesc(TTV_c(1:end,showperiod).*26.3163,'XData',time,clims);
title('control');
nexttile
imagesc(TTV_p(1:end,showperiod).*26.3163,'XData',time,clims);
title('NTS - PVH');
nexttile
imagesc(TTV_v(1:end,showperiod).*26.3163,'XData',time,clims);
title('NTS - VLM');


cb = colorbar;
cb.Layout.Tile = 'east';
title(tcl,' MAP deltaBP single trial, bin=1s, stimulation 0-60s ', fontsize=12);
%% Heatmap plot total volumne_averaged
 tt=40;
 showperiod = 51:90

bintime = 1;
nstim = 10; % time before stimulation
nbin = nstim/bintime;
%st =(-nbin+0.5)*bintime;
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
imagesc(mTTV_c(1:end,showperiod).*26.3163,'XData',time,clims);
title('Control');
nexttile
imagesc(mTTV_p(1:end,showperiod).*26.3163,'XData',time,clims);
title('NTS - PVH');
nexttile
imagesc(mTTV_v(1:end,showperiod).*26.3163,'XData',time,clims);
title('NTS - VLM');
cb = colorbar;
cb.Layout.Tile = 'east';
title(tcl,' delta BP mouse average 3 trials, bin=5s, stimulation 0-30s ', fontsize=12);


%% Heatmap plot total volumne_averaged _ RANKORDERED/Sort
 % tt = 135;
 tt=40;
 showperiod = 51:90

bintime = 1;
nstim = 10; % time before stimulation
nbin = nstim/bintime;
%st =(-nbin+0.5)*bintime;
st =(-nbin+0.5)*bintime;
et =tt-(nbin+0.5)*bintime;
% st =(-nbin+0.5)*bintime;
% et =tt-(nbin+0.5)*bintime;
time = st:bintime:et;
figure;
clims =[-10,10]

dataC=mTTV_c(1:end,:)
% Rank rows by their mean values (descending)
rowMeans = mean(dataC(:,61:85), 2);
[~, rowOrder] = sort(rowMeans, 'descend');
data_sorted_C = dataC(rowOrder, :).*26.3163;

dataP=mTTV_p(1:end,:)
% Rank rows by their mean values (descending)
rowMeans = mean(dataP(:,61:85), 2);
[~, rowOrder] = sort(rowMeans, 'descend');
data_sorted_P = dataP(rowOrder, :).*26.3163;

dataV=mTTV_v(1:end,:)
% Rank rows by their mean values (descending)
rowMeans = mean(dataV(:,61:85), 2);
[~, rowOrder] = sort(rowMeans, 'descend');
data_sorted_V = dataV(rowOrder, :).*26.3163;

colormap('jet')
tcl=tiledlayout(3,1);
nexttile
imagesc(data_sorted_C(:,showperiod ),'XData',time,clims);
title('Control');
nexttile
imagesc(data_sorted_P(:,showperiod ),'XData',time,clims);
title('NTS - PVH');
nexttile
imagesc(data_sorted_V(:,showperiod),'XData',time,clims);
title('NTS - VLM');

cb = colorbar;
cb.Layout.Tile = 'east';
title(tcl,' all 4 group delta BP mouse average 3 trials, bin=5s, stimulation 0-30s ', fontsize=12);
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
    lowbound = -5
    topbound = 2
    exchange_list='0123456789ABCDEF#'
    %colors ={'#0072BD','#D95319','#EDB120',	'#7E2F8E',	'#77AC30','#4DBEEE','#A2142F','#FF69B4'}
    %Colors for small volume
    % colors = { '#E274A9', '#13AF68','#0F80FF','#424242'};  % line colors %"#E274A9" "#26A7E1'#FFE009'
    colors = { '#424242','#0000FF','#FF0000'};  % line colors %"#E274A9" "#26A7E1'#FFE009'

    x=time
    h = figure('visible','on');
    hA = axes(h);

    mm=[]
    error=[]
    var = data_sorted_C(:,showperiod );
    mm(:,1) = mean(var,'omitnan')
    error(:,1)= std(var,'omitnan')/sqrt(size(var,1));
    var = data_sorted_P(:,showperiod );
    mm(:,2) = mean(var,'omitnan')
    error(:,2)= std(var,'omitnan')/sqrt(size(var,1));
    var = data_sorted_V(:,showperiod );
    mm(:,3) = mean(var,'omitnan')
    error(:,3)= std(var,'omitnan')/sqrt(size(var,1));

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
   
    group = {'Control','NTS-PVH','NTS-VLM'}
    lgd=legend(group,'Location','NorthEastOutside','Box','off')
    title(lgd,'Group')
    
    %fill([0,0,stim,stim],[lowbound,topbound,topbound,lowbound],'yellow','FaceAlpha',0.2,'linestyle','none');



%% BP-Select 2*2*3 representative trace
j=0
o=0
l=0
m=0
downsample=100 %downsample times
time = -29.999:0.001* downsample:60; 
snn= (30+0.001)*(sampleRate);
enn = 120*(sampleRate);
showwindow = int32([snn:1:enn]);
figure
%colors = { '#E274A9','#E274A9', '#E274A9','#E274A9','#E274A9','#E274A9','#13AF68','#13AF68','#0F80FF','#0F80FF','#424242','#424242'};  % line colors %"#E274A9" "#26A7E1'#FFE009'
colors =  {"#424242","#424242","#424242","#424242","#0000FF","#0000FF","#0000FF","#0000FF","#FF0000","#FF0000","#FF0000","#FF0000","#FF0000"};
ii=0
tracelabel={}
%ii = [2 18 22 29]
%jj = [1 2 2 1]
ii = [38 39 71 70 63 64 19 21 28 29 32 33]

cd = data_extract{2,3};
ci = trialinfo_extract{2,3};
for i = 1:length(ii)
    mm = ii(i);
    current=cd(showwindow,mm);
    y0 = smoothdata(current,'movmean',10);
    y = decimate(y0, downsample);
    plot(time,y-l-max(y),'Color',colors{i});
    l = l+(max(y)-min(y))*1.1;
    o = o+1;
    tracelabel{o} = [ci{mm,4},'-',ci{mm,3},'-',ci{mm,2},'-T',num2str(ci{mm,5}),'-idx',num2str(mm)];
    hold on
           
end
legend(tracelabel);
xregion(0, 60);

plot([-5 -5], [2 2.076],'k-',[-5 0],[2 2],'k-','linewidth',3); % 3point for the bar
text(-4.5,1.5,'5s','horiz','center','vert','top'); 
text(-5.5,4.5,'2cmH2O','horiz','right','vert','top')
% plot([-1 -1 0],[2.1 2 2],'k-','linewidth',3); % 3point for the bar
% text(-0.5,1.5,'1s','horiz','center','vert','top'); 
% text(80,.85,'0.1 y units ','horiz','right','vert','middle');
fname ='BP-1min-all-example trace'
title(fname, fontsize=12);

print(gcf,'-vector','-dsvg',[fname,'.svg']) % svg
%print(gcf,'-vector','-dpdf',[fname,'.pdf']) % pdf
