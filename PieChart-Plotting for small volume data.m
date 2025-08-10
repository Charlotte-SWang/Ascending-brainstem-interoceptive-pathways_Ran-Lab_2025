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
cc = single_z_amp_tv


info_current = trialinfo_extract{2,3};
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
  %% Heatmap plot BR
    tt = 15;
    bintime = 1;
    stim = 5;
    nbin = stim/bintime;
    st =(-nbin+1)*bintime;
    et =tt-(nbin)*bintime;
    time = st:bintime:et;
    showperiod = [1:15]
clims = [-5 5]
colormap('jet')
tcl=tiledlayout(3,1);
nexttile
imagesc(TTV_c(1:end,showperiod),'XData',time,clims);
title('control');
nexttile
imagesc(TTV_p(1:end,showperiod),'XData',time,clims);
title('NTS - PVH');
nexttile
imagesc(TTV_v(1:end,showperiod),'XData',time,clims);
title('NTS - VLM');


cb = colorbar;
cb.Layout.Tile = 'east';
title(tcl,' BR zscore single trial, bin=5s, stimulation 0-5s ', fontsize=12);

%% Heatmap plot BR_averaged
    tt = 15;
    bintime = 1;
    stim = 5;
    nbin = stim/bintime;
    st =(-nbin+0.5)*bintime;
    et =tt-(nbin+0.5)*bintime;
    time = st:bintime:et;
    showperiod = [1:15]
clims = [-5 5]

colormap('jet')
tcl=tiledlayout(3,1);
nexttile
imagesc(mTTV_c(1:end,showperiod),'XData',time,clims);
title('Control');
nexttile
imagesc(mTTV_p(1:end,showperiod),'XData',time,clims);
title('NTS - PVH');
nexttile
imagesc(mTTV_v(1:end,showperiod),'XData',time,clims);
title('NTS - VLM');
cb = colorbar;
cb.Layout.Tile = 'east';
title(tcl,' delta BP mouse average 3 trials, bin=5s, stimulation 0-25s ', fontsize=12);

%% Plot Breathing (bin =1s=)
    tt = 15;
    bintime = 1;
    stim = 5;
    nbin = stim/bintime;
    st =(-nbin)*bintime;
    et =tt-(nbin)*bintime;
    time = st:bintime:et;
    showperiod = [1:15]
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
   
    group = {'Anterior','intermediate','posterior','Control'}
    lgd=legend(group,'Location','NorthEastOutside','Box','off')
    title(lgd,'Group')
    
    %fill([0,0,stim,stim],[lowbound,topbound,topbound,lowbound],'yellow','FaceAlpha',0.2,'linestyle','none');

%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
 %% plot BR value
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
cc = BR_CI


info_current = trialinfo_extract{2,3};
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
  %% Heatmap plot BR
    tt = 15;
    bintime = 1;
    stim = 5;
    nbin = stim/bintime;
    st =(-nbin+0.5)*bintime;
    et =tt-(nbin+0.5)*bintime;
    time = st:bintime:et;
    showperiod = [1:15]
clims = [-0.5 0.5]
colormap('jet')
tcl=tiledlayout(3,1);
nexttile
imagesc(TTV_c(1:end,showperiod),'XData',time,clims);
title('control');
nexttile
imagesc(TTV_p(1:end,showperiod),'XData',time,clims);
title('NTS - PVH');
nexttile
imagesc(TTV_v(1:end,showperiod),'XData',time,clims);
title('NTS - VLM');


cb = colorbar;
cb.Layout.Tile = 'east';
title(tcl,' BR zscore single trial, bin=5s, stimulation 0-5s ', fontsize=12);

%% Heatmap plot BR_averaged
    tt = 15;
    bintime = 1;
    stim = 5;
    nbin = stim/bintime;
    st =(-nbin+0.5)*bintime;
    et =tt-(nbin+0.5)*bintime;
    time = st:bintime:et;
    showperiod = [1:15]
clims = [-5 5]

colormap('jet')
tcl=tiledlayout(3,1);
nexttile
imagesc(mTTV_c(1:end,showperiod),'XData',time,clims);
title('Control');
nexttile
imagesc(mTTV_p(1:end,showperiod),'XData',time,clims);
title('NTS - PVH');
nexttile
imagesc(mTTV_v(1:end,showperiod),'XData',time,clims);
title('NTS - VLM');
cb = colorbar;
cb.Layout.Tile = 'east';
title(tcl,' delta BP mouse average 3 trials, bin=5s, stimulation 0-25s ', fontsize=12);

%% Plot Breathing (bin =1s=)
    tt = 15;
    bintime = 1;
    stim = 5;
    nbin = stim/bintime;
    st =(-nbin)*bintime;
    et =tt-(nbin)*bintime;
    time = st:bintime:et;
    showperiod = [1:15]
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
   
    group = {'Anterior','intermediate','posterior','Control'}
    lgd=legend(group,'Location','NorthEastOutside','Box','off')
    title(lgd,'Group')
    
    %fill([0,0,stim,stim],[lowbound,topbound,topbound,lowbound],'yellow','FaceAlpha',0.2,'linestyle','none');


%% Z-score analysis for Tidal Volume
  passz = 2;
  passbnumber = 3;
  passtrialnumber = 2;

  tv_z_conclusion = [];

  % info_current = trialinfo_extract{3,2};
  % [unique_vals, ~, idx] = unique([info_current{:,1}],'stable');
  % counts = accumarray(idx, 1);  % returns 5
  % 
  zz = single_z_amp_tv;
  nrep=3;
  for k = 1:size(info,1)
    %ccsubset =cc((sum(counts(1:k))-counts(k)+1):sum(counts(1:k)))
    %ccmat= cat(1,ccsubset{:})
    y_t_increase=0;
    y_t_decrease=0;
    for nn = 1:nrep
        kk = (k-1)*nrep +nn;
          % passpercent = 0.5;
          % passnumber=floor(length(zz{kk})*0.05);
          % if passnumber == 0
          %     passnumber = 1
          % end
        zz_max = maxk(zz{kk},passbnumber);
        zz_min = mink(zz{kk},passbnumber);
        if min(zz_max)>=passz
            y_t_increase = y_t_increase+1;
        end
        if abs(max(zz_min))>=passz & max(zz_min)<0
            y_t_decrease = y_t_decrease+1;
        end
    end
    if y_t_increase >= passtrialnumber
        tv_z_conclusion(k,1)=1;
    else
        tv_z_conclusion(k,1)=0;
    end
     if y_t_decrease >= passtrialnumber
        tv_z_conclusion(k,2)=1;
    else
        tv_z_conclusion(k,2)=0;
     end
  end