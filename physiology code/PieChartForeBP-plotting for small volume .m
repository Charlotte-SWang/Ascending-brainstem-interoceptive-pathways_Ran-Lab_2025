%% plot BP zscore
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
cc = MAP_zscore

info_current = trialinfo_extract{2,3};
[unique_vals, ~, idx] = unique([info_current{:,1}]);
counts = accumarray(idx, 1);  % returns 5
cc_zscore=[];

for k = 1:size(counts,1)
    ccsubset =cc((sum(counts(1:k))-counts(k)+1):sum(counts(1:k)));
    ccmat= cat(1,ccsubset{:});
    ccmat_w= cat(2,ccsubset{:});
    ccmat_window= mean(ccmat_w(61:70,:),1);
    if size(ccmat_window,2)==2
         ccmat_window= [ccmat_window, NaN];
    end
    cc_zscore(k,:)=mean(ccmat_window,1,"omitnan");

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
  %% Heatmap plot Blood pressure
    tt = 73;
    bintime = 1;
    nstim = 10; % time before stimulation
    nbin = nstim/bintime;
    %st =(-nbin+0.5)*bintime;
    st =(-nbin+0.5)*bintime;
    et =tt-(nbin+0.5)*bintime;
    time = st:bintime:et;
    showperiod = [51:123]
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
title(tcl,' BP zscore single trial, bin=5s, stimulation 0-5s ', fontsize=12);
