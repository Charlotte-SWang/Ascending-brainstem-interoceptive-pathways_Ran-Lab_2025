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
colors =  {"#424242","#424242","#424242","#424242","#0000FF","#0000FF","#0000FF","#0000FF","#FF0000","#FF0000","#FF0000","#FF0000","#FF0000"};
ii=0
tracelabel={}
%ii = [2 18 22 29]
%jj = [1 2 2 1]
ii = [52 53 82 83 14 15 49 50 47 48 89 90]
cd = data_extract{3,1};
ci = trialinfo_extract{3,1};
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
fname ='HR-5s-example trace'
title(fname, fontsize=12);

print(gcf,'-vector','-dsvg',[fname,'.svg']) % svg
%print(gcf,'-vector','-dpdf',[fname,'.pdf']) % pdf