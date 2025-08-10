%% plot GP Total
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
cc = GP;


info_current = trialinfo_extract{1,4};
[unique_vals, ~, idx] = unique([info_current{:,1}],'stable');
counts = accumarray(idx, 1);  % returns 5

tt = 615;
bintime = 0.001;
stim = 300;
nbin = stim/bintime;
st =(-nbin+1)*bintime;
et =tt-(nbin)*bintime;
time = st:bintime:et;
%baseline = [226001:300000];
baseline = [1:300000];
windowSize = 1000;             % number of points in the moving window

for k = 1:size(counts,1)
    %ccsubset =cc(k,:);
    %ccmat= cat(1,ccsubset{:});
    ccsubset = {};
    for nn = 1:counts(k)
        if k==1
            kk = nn;
        else
            kk = sum(counts(1:k-1)) +nn;
        end
        % %%%% AUC
        % Fs = 1000;        % Sampling frequency (samples per second)
        % bin_sec = 15;     % Bin size in seconds
        % bin_size = Fs * bin_sec;   % Bin size in samples
        % 
        % t = (0:length(cc{k,nn})-1) / Fs;  % Time vector in seconds
        % n_bins = floor(length(cc{k,nn}) / bin_size);
        % AUCs = zeros(1, n_bins);
        % 
        % for i = 1:n_bins
        %     idx_start = (i-1)*bin_size + 1;
        %     idx_end = i*bin_size;
        %     AUCs(i) = trapz(t(idx_start:idx_end), cc{k,nn}(idx_start:idx_end));
        % end
        % cc{k,nn}=AUCs;

        % % % %%%%%% %FFT -Start
        % % % Fs = 1000;                 % Sampling frequency (Hz)
        % % % t = 0:1/Fs:300-1/Fs;         % Time vector (1 second)
        % % % x1 =  cc{k,nn}(1:300000);        % Example signal: 50 Hz sine wave
        % % % x2 =  cc{k,nn}(300001:600000);
        % % % x=x1;
        % % % Y=[]
        % % % Y = fft(x);                % Compute FFT
        % % % N = length(Y);             % Number of points
        % % % f = Fs*(0:(N/2))/N;        % Frequency vector (one-sided)
        % % % % Compute magnitude (one-sided)
        % % % P = abs(Y/N);              % Normalize
        % % % P1 = P(1:N/2+1);           % Take positive frequencies
        % % % P1(2:end-1) = 2*P1(2:end-1); % Double non-DC, non-Nyquist components
        % % % 
        % % % x=x2;
        % % % Y = fft(x);                % Compute FFT
        % % % N = length(Y);             % Number of points
        % % % f2 = Fs*(0:(N/2))/N;        % Frequency vector (one-sided)
        % % % % Compute magnitude (one-sided)
        % % % P = abs(Y/N);              % Normalize
        % % % P2 = P(1:N/2+1);           % Take positive frequencies
        % % % P2(2:end-1) = 2*P2(2:end-1); % Double non-DC, non-Nyquist components
        % % % 
        % % % % Plot frequency spectrum
        % % % figure;
        % % % plot(f, P1);
        % % % hold on;
        % % % plot(f2,P2);
        % % % xlabel('Frequency (Hz)');
        % % % ylabel('|P(f)|');
        % % % xlim([0,0.5])
        % % % title('Single-Sided Amplitude Spectrum');
        % % % %%%%%% %FFT -End
        

        %%% Running average
        cc{k,nn} = movmean(cc{k,nn}, windowSize);
        cc{k,nn} =cc{k,nn} -min(cc{k,nn} )
        %%%%%Normalize
        %baseline = 1:20;
        cc_z= mean(cc{k,nn}(baseline));
        cc_sd = std(cc{k,nn}(baseline));
        cc_nn = (cc{k,nn}-cc_z)./cc_sd;
        cc_nn = (cc{k,nn}-cc_z);
        % cc_min = mean(cc{k,nn}(baseline));
        % cc_nn = (cc{k,nn})-cc_min;

        if ismember('C',info_current{kk,4})
            %TTV_c(j,:) = cc{k,nn};
            TTV_c(j,:) = cc_nn;
            j = j+1;
        else
            if ismember('PVH',info_current{kk,3})
                %TTV_p(o,:) = cc{k,nn};
                TTV_p(o,:) = cc_nn;
                o = o+1;
            elseif ismember('VLM',info_current{kk,3})
                % TTV_v(l,:) = cc{k,nn};
                TTV_v(l,:) = cc_nn;
                l = l+1;
            end
        end
        ccsubset = [ccsubset; cc_nn];
 
    end
    
    catidx=2;
    ccmat= cat(catidx,ccsubset{:});

    if ismember('C',info_current{kk,4})
        mTTV_c(jj,:) = sum(ccmat,catidx,"omitnan")/counts(k);
        jj=jj+1;
    else
        if ismember('PVH',info_current{kk,3})
            mTTV_p(oo,:) =  sum(ccmat,catidx,"omitnan")/counts(k);
            oo=oo+1;
        elseif ismember('VLM',info_current{kk,3})
            mTTV_v(ll,:) = sum(ccmat,catidx,"omitnan")/counts(k);
            ll=ll+1;
        end
    end
end
  %% Heatmap plot 
    tt = 615;
    bintime = 0.001;
    stim = 300;
    nbin = stim/bintime;
    st =(-nbin+1)*bintime;
    et =tt-(nbin)*bintime;
    time = st:bintime:et;
    showperiod = [1:1:615000];
    %showperiod = [1:1:41];
clims = [-2 2];
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
title(tcl,' GP total substract zscore single trial, bin=15s, stimulation 0-300s ', fontsize=12);

%% Heatmap plot total volumne_averaged _ RANKORDERED/Sort
rowOrder = [];
baseline = [300001:600000];%
%baseline =[21:40]
clims = [-5 5];

% [300001:420000];
dataC=mTTV_c(1:end,:);
%clear sort
% Rank rows by their mean values (descending)
rowMeans = mean(dataC(:,baseline), 2);

[~, rowOrder] = sort(rowMeans, 'descend');
data_sorted_C = dataC(rowOrder, :);

dataP=mTTV_p(1:end,:);
% Rank rows by their mean values (descending)
rowMeans = mean(dataP(:,baseline), 2);
[~, rowOrder] = sort(rowMeans, 'descend');
data_sorted_P = dataP(rowOrder, :);

dataV=mTTV_v(1:end,:);
% Rank rows by their mean values (descending)
rowMeans = mean(dataV(:,baseline), 2);
[~, rowOrder] = sort(rowMeans, 'descend');
data_sorted_V = dataV(rowOrder, :);

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
title(tcl,' zscore GP tonic mouse average 3 trials, bin=15s, stimulation 0-300s ', fontsize=12);

%% plot GP
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
cc = GP_Tonic;


info_current = trialinfo_extract{1,4};
[unique_vals, ~, idx] = unique([info_current{:,1}],'stable');
counts = accumarray(idx, 1);  % returns 5

tt = 615;
bintime = 0.001;
stim = 300;
nbin = stim/bintime;
st =(-nbin+1)*bintime;
et =tt-(nbin)*bintime;
time = st:bintime:et;
%baseline = [226001:300000];
baseline = [1:300000];
windowSize = 15000;             % number of points in the moving window

for k = 1:size(counts,1)
    %ccsubset =cc(k,:);
    %ccmat= cat(1,ccsubset{:});
    ccsubset = {};
    for nn = 1:counts(k)
        if k==1
            kk = nn;
        else
            kk = sum(counts(1:k-1)) +nn;
        end
        % %%%% AUC
        % Fs = 1000;        % Sampling frequency (samples per second)
        % bin_sec = 15;     % Bin size in seconds
        % bin_size = Fs * bin_sec;   % Bin size in samples
        % 
        % t = (0:length(cc{k,nn})-1) / Fs;  % Time vector in seconds
        % n_bins = floor(length(cc{k,nn}) / bin_size);
        % AUCs = zeros(1, n_bins);
        % 
        % for i = 1:n_bins
        %     idx_start = (i-1)*bin_size + 1;
        %     idx_end = i*bin_size;
        %     AUCs(i) = trapz(t(idx_start:idx_end), cc{k,nn}(idx_start:idx_end));
        % end
        % cc{k,nn}=AUCs;

        % % % %%%%%% %FFT -Start
        % % % Fs = 1000;                 % Sampling frequency (Hz)
        % % % t = 0:1/Fs:300-1/Fs;         % Time vector (1 second)
        % % % x1 =  cc{k,nn}(1:300000);        % Example signal: 50 Hz sine wave
        % % % x2 =  cc{k,nn}(300001:600000);
        % % % x=x1;
        % % % Y=[]
        % % % Y = fft(x);                % Compute FFT
        % % % N = length(Y);             % Number of points
        % % % f = Fs*(0:(N/2))/N;        % Frequency vector (one-sided)
        % % % % Compute magnitude (one-sided)
        % % % P = abs(Y/N);              % Normalize
        % % % P1 = P(1:N/2+1);           % Take positive frequencies
        % % % P1(2:end-1) = 2*P1(2:end-1); % Double non-DC, non-Nyquist components
        % % % 
        % % % x=x2;
        % % % Y = fft(x);                % Compute FFT
        % % % N = length(Y);             % Number of points
        % % % f2 = Fs*(0:(N/2))/N;        % Frequency vector (one-sided)
        % % % % Compute magnitude (one-sided)
        % % % P = abs(Y/N);              % Normalize
        % % % P2 = P(1:N/2+1);           % Take positive frequencies
        % % % P2(2:end-1) = 2*P2(2:end-1); % Double non-DC, non-Nyquist components
        % % % 
        % % % % Plot frequency spectrum
        % % % figure;
        % % % plot(f, P1);
        % % % hold on;
        % % % plot(f2,P2);
        % % % xlabel('Frequency (Hz)');
        % % % ylabel('|P(f)|');
        % % % xlim([0,0.5])
        % % % title('Single-Sided Amplitude Spectrum');
        % % % %%%%%% %FFT -End
        

        %%% Running average
        cc{k,nn} = movmean(cc{k,nn}, windowSize);

        %%%%%Normalize
        %baseline = 1:20;
        cc_z= mean(cc{k,nn}(baseline));
        cc_sd = std(cc{k,nn}(baseline));
        cc_nn = (cc{k,nn}-cc_z)./cc_sd;

        % cc_min = mean(cc{k,nn}(baseline));
        % cc_nn = (cc{k,nn})-cc_min;

        if ismember('C',info_current{kk,4})
            %TTV_c(j,:) = cc{k,nn};
            TTV_c(j,:) = cc_nn;
            j = j+1;
        else
            if ismember('PVH',info_current{kk,3})
                %TTV_p(o,:) = cc{k,nn};
                TTV_p(o,:) = cc_nn;
                o = o+1;
            elseif ismember('VLM',info_current{kk,3})
                % TTV_v(l,:) = cc{k,nn};
                TTV_v(l,:) = cc_nn;
                l = l+1;
            end
        end
        ccsubset = [ccsubset; cc_nn];
 
    end
    
    catidx=2;
    ccmat= cat(catidx,ccsubset{:});

    if ismember('C',info_current{kk,4})
        mTTV_c(jj,:) = sum(ccmat,catidx,"omitnan")/counts(k);
        jj=jj+1;
    else
        if ismember('PVH',info_current{kk,3})
            mTTV_p(oo,:) =  sum(ccmat,catidx,"omitnan")/counts(k);
            oo=oo+1;
        elseif ismember('VLM',info_current{kk,3})
            mTTV_v(ll,:) = sum(ccmat,catidx,"omitnan")/counts(k);
            ll=ll+1;
        end
    end
end
  %% Heatmap plot 
    tt = 615;
    bintime = 0.001;
    stim = 300;
    nbin = stim/bintime;
    st =(-nbin+1)*bintime;
    et =tt-(nbin)*bintime;
    time = st:bintime:et;
    showperiod = [1:1:615000];
    %showperiod = [1:1:41];
clims = [-5 5];
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
title(tcl,' GP tonic zscore single trial, bin=15s, stimulation 0-300s ', fontsize=12);

%% Heatmap plot total volumne_averaged _ RANKORDERED/Sort
rowOrder = [];
baseline = [300001:600000];%
%baseline =[21:40]
clims = [-5 5];

% [300001:420000];
dataC=mTTV_c(1:end,:);
%clear sort
% Rank rows by their mean values (descending)
rowMeans = mean(dataC(:,baseline), 2);

[~, rowOrder] = sort(rowMeans, 'descend');
data_sorted_C = dataC(rowOrder, :);

dataP=mTTV_p(1:end,:);
% Rank rows by their mean values (descending)
rowMeans = mean(dataP(:,baseline), 2);
[~, rowOrder] = sort(rowMeans, 'descend');
data_sorted_P = dataP(rowOrder, :);

dataV=mTTV_v(1:end,:);
% Rank rows by their mean values (descending)
rowMeans = mean(dataV(:,baseline), 2);
[~, rowOrder] = sort(rowMeans, 'descend');
data_sorted_V = dataV(rowOrder, :);

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
title(tcl,' zscore GP tonic mouse average 3 trials, bin=15s, stimulation 0-300s ', fontsize=12);

%% Plot GP (bin =15s=)
    tt = 615;
    bintime = 0.001;
    stim =300;
    nbin = stim/bintime;
    st =(-nbin+1)*bintime;
    et =tt-(nbin)*bintime;
    time = st:bintime:et;
    %showperiod = [1:41];
    showperiod = [1:1:615000];
    % ylim([lowbound,topbound])
    % xlim([(-nbin+1)*bintime (tt/bin-nbin)*bintime])
    % set(gca,'FontSize', 15);
    % hold on
    % 
    % fill([0,0,120,120],[lowbound,topbound,topbound,lowbound],'yellow','FaceAlpha',0.2,'linestyle','none');
    % 
    lowbound = -10;
    topbound = 10;
    exchange_list='0123456789ABCDEF#';
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
    mm(:,1) = mean(var,'omitnan');
    error(:,1)= std(var,'omitnan')/sqrt(size(var,1));
    var = data_sorted_P(:,showperiod );
    mm(:,2) = mean(var,'omitnan');
    error(:,2)= std(var,'omitnan')/sqrt(size(var,1));
    var = data_sorted_V(:,showperiod );
    mm(:,3) = mean(var,'omitnan');
    error(:,3)= std(var,'omitnan')/sqrt(size(var,1));

    count=1
    p=[]
    for i=1:size(mm,2)
        y1=mm(:,i)+error(:,i);
        y2=mm(:,i)-error(:,i);
        x2=[x,fliplr(x)];
        inBetween = [y1',fliplr(y2')];
        string=(colors{count});
        num=zeros(1,3);
        for r=1:3
            temp_Coe1=find(exchange_list==string(r*2))-1;
            temp_Coe2=find(exchange_list==string(r*2+1))-1;
            num(r)=(16*temp_Coe1+temp_Coe2)/255;
        end
        fi=fill(x2,inBetween,num,'HandleVisibility','off');
        set(fi,'edgealpha',0,'facealpha',0.3);
        hold on
        plot(x,mm(:,i), 'LineWidth',1, 'Color', colors{count});

        count=count+1;
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
   
    group = {'Control','NTS-PVH','NTS-VLM'};
    lgd=legend(group,'Location','NorthEastOutside','Box','off');
    title(lgd,'Group');
    
    
    

    
    
    
    
    
    
%% plot GP Phasic
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
cc = GP_Phasic;


info_current = trialinfo_extract{1,4};
[unique_vals, ~, idx] = unique([info_current{:,1}],'stable');
counts = accumarray(idx, 1);  % returns 5

tt = 615;
bintime = 0.001;
stim = 300;
nbin = stim/bintime;
st =(-nbin+1)*bintime;
et =tt-(nbin)*bintime;
time = st:bintime:et;
%baseline = [226001:300000];
baseline = [1:300000];
windowSize = 15000;             % number of points in the moving window

for k = 1:size(counts,1)
    %ccsubset =cc(k,:);
    %ccmat= cat(1,ccsubset{:});
    ccsubset = {};
    for nn = 1:counts(k)
        if k==1
            kk = nn;
        else
            kk = sum(counts(1:k-1)) +nn;
        end
        %%%% AUC
        Fs = 1000;        % Sampling frequency (samples per second)
        bin_sec = 15;     % Bin size in seconds
        bin_size = Fs * bin_sec;   % Bin size in samples

        t = (0:length(cc{k,nn})-1) / Fs;  % Time vector in seconds
        n_bins = floor(length(cc{k,nn}) / bin_size);
        AUCs = zeros(1, n_bins);

        for i = 1:n_bins
            idx_start = (i-1)*bin_size + 1;
            idx_end = i*bin_size;
            AUCs(i) = trapz(t(idx_start:idx_end), cc{k,nn}(idx_start:idx_end));
        end
        cc{k,nn}=AUCs;

        % % % %%%%%% %FFT -Start
        % % % Fs = 1000;                 % Sampling frequency (Hz)
        % % % t = 0:1/Fs:300-1/Fs;         % Time vector (1 second)
        % % % x1 =  cc{k,nn}(1:300000);        % Example signal: 50 Hz sine wave
        % % % x2 =  cc{k,nn}(300001:600000);
        % % % x=x1;
        % % % Y=[]
        % % % Y = fft(x);                % Compute FFT
        % % % N = length(Y);             % Number of points
        % % % f = Fs*(0:(N/2))/N;        % Frequency vector (one-sided)
        % % % % Compute magnitude (one-sided)
        % % % P = abs(Y/N);              % Normalize
        % % % P1 = P(1:N/2+1);           % Take positive frequencies
        % % % P1(2:end-1) = 2*P1(2:end-1); % Double non-DC, non-Nyquist components
        % % % 
        % % % x=x2;
        % % % Y = fft(x);                % Compute FFT
        % % % N = length(Y);             % Number of points
        % % % f2 = Fs*(0:(N/2))/N;        % Frequency vector (one-sided)
        % % % % Compute magnitude (one-sided)
        % % % P = abs(Y/N);              % Normalize
        % % % P2 = P(1:N/2+1);           % Take positive frequencies
        % % % P2(2:end-1) = 2*P2(2:end-1); % Double non-DC, non-Nyquist components
        % % % 
        % % % % Plot frequency spectrum
        % % % figure;
        % % % plot(f, P1);
        % % % hold on;
        % % % plot(f2,P2);
        % % % xlabel('Frequency (Hz)');
        % % % ylabel('|P(f)|');
        % % % xlim([0,0.5])
        % % % title('Single-Sided Amplitude Spectrum');
        % % % %%%%%% %FFT -End
        

        % %%% Running average
        % cc{k,nn} = movmean(cc{k,nn}, windowSize);

        %%%%%Normalize
        baseline = 1:20;
        cc_z= mean(cc{k,nn}(baseline));
        cc_sd = std(cc{k,nn}(baseline));
        cc_nn = (cc{k,nn}-cc_z)./cc_sd;

        % cc_min = mean(cc{k,nn}(baseline));
        % cc_nn = (cc{k,nn})-cc_min;

        if ismember('C',info_current{kk,4})
            %TTV_c(j,:) = cc{k,nn};
            TTV_c(j,:) = cc_nn;
            j = j+1;
        else
            if ismember('PVH',info_current{kk,3})
                %TTV_p(o,:) = cc{k,nn};
                TTV_p(o,:) = cc_nn;
                o = o+1;
            elseif ismember('VLM',info_current{kk,3})
                % TTV_v(l,:) = cc{k,nn};
                TTV_v(l,:) = cc_nn;
                l = l+1;
            end
        end
        ccsubset = [ccsubset; cc_nn];
 
    end
    
    catidx=1;
    ccmat= cat(catidx,ccsubset{:});

    if ismember('C',info_current{kk,4})
        mTTV_c(jj,:) = sum(ccmat,catidx,"omitnan")/counts(k);
        jj=jj+1;
    else
        if ismember('PVH',info_current{kk,3})
            mTTV_p(oo,:) =  sum(ccmat,catidx,"omitnan")/counts(k);
            oo=oo+1;
        elseif ismember('VLM',info_current{kk,3})
            mTTV_v(ll,:) = sum(ccmat,catidx,"omitnan")/counts(k);
            ll=ll+1;
        end
    end
end
  %% Heatmap plot 
    tt = 615;
    bintime = 15;
    stim = 300;
    nbin = stim/bintime;
    st =(-nbin+1)*bintime;
    et =tt-(nbin)*bintime;
    time = st:bintime:et;
    %showperiod = [1:1:615000];
    showperiod = [1:1:41];
clims = [-5 5];
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
title(tcl,' GP phasic zscore single trial, bin=15s, stimulation 0-300s ', fontsize=12);

%% Heatmap plot total volumne_averaged _ RANKORDERED/Sort
rowOrder = [];
baseline = [300001:600000];%
baseline =[21:40]
clims = [-5 5];

% [300001:420000];
dataC=mTTV_c(1:end,:);
%clear sort
% Rank rows by their mean values (descending)
rowMeans = mean(dataC(:,baseline), 2);

[~, rowOrder] = sort(rowMeans, 'descend');
data_sorted_C = dataC(rowOrder, :);

dataP=mTTV_p(1:end,:);
% Rank rows by their mean values (descending)
rowMeans = mean(dataP(:,baseline), 2);
[~, rowOrder] = sort(rowMeans, 'descend');
data_sorted_P = dataP(rowOrder, :);

dataV=mTTV_v(1:end,:);
% Rank rows by their mean values (descending)
rowMeans = mean(dataV(:,baseline), 2);
[~, rowOrder] = sort(rowMeans, 'descend');
data_sorted_V = dataV(rowOrder, :);

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
title(tcl,' zscore GP phasic mouse average 3 trials, bin=15s, stimulation 0-300s ', fontsize=12);

%% Plot GP (bin =15s=)
    tt = 615;
    bintime = 15;
    stim =300;
    nbin = stim/bintime;
    st =(-nbin+1)*bintime;
    et =tt-(nbin)*bintime;
    time = st:bintime:et;
    showperiod = [1:41];
    %showperiod = [1:1:615000];
    % ylim([lowbound,topbound])
    % xlim([(-nbin+1)*bintime (tt/bin-nbin)*bintime])
    % set(gca,'FontSize', 15);
    % hold on
    % 
    % fill([0,0,120,120],[lowbound,topbound,topbound,lowbound],'yellow','FaceAlpha',0.2,'linestyle','none');
    % 
    lowbound = -2;
    topbound = 2;
    exchange_list='0123456789ABCDEF#';
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
    mm(:,1) = mean(var,'omitnan');
    error(:,1)= std(var,'omitnan')/sqrt(size(var,1));
    var = data_sorted_P(:,showperiod );
    mm(:,2) = mean(var,'omitnan');
    error(:,2)= std(var,'omitnan')/sqrt(size(var,1));
    var = data_sorted_V(:,showperiod );
    mm(:,3) = mean(var,'omitnan');
    error(:,3)= std(var,'omitnan')/sqrt(size(var,1));

    count=1
    p=[]
    for i=1:size(mm,2)
        y1=mm(:,i)+error(:,i);
        y2=mm(:,i)-error(:,i);
        x2=[x,fliplr(x)];
        inBetween = [y1',fliplr(y2')];
        string=(colors{count});
        num=zeros(1,3);
        for r=1:3
            temp_Coe1=find(exchange_list==string(r*2))-1;
            temp_Coe2=find(exchange_list==string(r*2+1))-1;
            num(r)=(16*temp_Coe1+temp_Coe2)/255;
        end
        fi=fill(x2,inBetween,num,'HandleVisibility','off');
        set(fi,'edgealpha',0,'facealpha',0.3);
        hold on
        plot(x,mm(:,i), 'LineWidth',1, 'Color', colors{count});

        count=count+1;
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
   
    group = {'Control','NTS-PVH','NTS-VLM'};
    lgd=legend(group,'Location','NorthEastOutside','Box','off');
    title(lgd,'Group');


    
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



%% GP-Select 2*2*3 representative trace
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
colors =  {"#424242","#424242","#424242","#424242","#424242","#424242","#0000FF","#0000FF","#0000FF","#0000FF","#0000FF","#0000FF","#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#FF0000"};
ii=0
tracelabel={}
%ii = [2 18 22 29]
%jj = [1 2 2 1]
ii = [93 94 105 106 41 42 23 26 20 21 11 12 56 57 83 85 87 88]

cd = data_extract{1,4};
ci = trialinfo_extract{1,4};
for i = 1:length(ii)
    mm = ii(i);
    plot(time,cd(showwindow,mm)-o*10,'Color',colors{i});
    o = o+1;
    tracelabel{o} = [ci{mm,4},'-',ci{mm,3},'-',ci{mm,2},'-T',num2str(ci{mm,5}),'-idx',num2str(mm)];
    hold on
           

end
legend(tracelabel);
xregion(0, 300);

plot([-5 -5], [2 3],'k-',[-30 0],[2 2],'k-','linewidth',3); % 3point for the bar
text(-4.5,1.5,'30s','horiz','center','vert','top'); 
text(-5.5,4.5,'1cmH2O','horiz','right','vert','top')
% plot([-1 -1 0],[2.1 2 2],'k-','linewidth',3); % 3point for the bar
% text(-0.5,1.5,'1s','horiz','center','vert','top'); 
% text(80,.85,'0.1 y units ','horiz','right','vert','middle');
fname ='GP-5min-example trace'
title(fname, fontsize=12);

print(gcf,'-vector','-dsvg',[fname,'.svg']) % svg
%print(gcf,'-vector','-dpdf',[fname,'.pdf']) % pdf
