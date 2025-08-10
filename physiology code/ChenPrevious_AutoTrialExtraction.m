%% Extract trial information
% Shiqi Wang 20250114
 
clear all


%% Input

%%% Attention !!!%%%
SampleRate = 1; % default unit: ms --check 'isi' & 'isi_units'

% Read instructions for protocol extraction
[ProtocolName,ProtocolPath]=uigetfile('*.csv');
cd(ProtocolPath);
protocol_set = readmatrix(ProtocolName);

%%% Need to check
FixInfo = 6 % Channel is from 7:end

%Read summary for all the files
[FileSummary,FileSummaryPath]=uigetfile('*.xlsx');
cd(FileSummaryPath);
individual_protocol = readtable(FileSummary,ReadRowNames=true);

% Get all data
[pathname1] = uigetdir( 'Select the physiology data');
cd(pathname1);
filename_sorted=dir('*-*-*-*.mat');
n=length(filename_sorted);%file numbers
data_extract = {};
trialinfo_extract={};
trial_order=zeros(3,4);
for i = 1:n
    % Prepare information
    [~,token,~] = fileparts(filename_sorted(i).name);
    token_split = split(token,'-');

    order = str2num(token_split{1});


    info{order,1} = order;
    info{order,2} = token_split{2};
    info{order,3}= token_split{3};
    info{order,4} = token_split{4};

    % Read data file
    [~, rowNum] = ismember(filename_sorted(i).name, individual_protocol.Properties.RowNames);
    extract_guide = individual_protocol{rowNum,:};
    %%m = max(extract_guide);
    
    % Read data
    fileDirName = [pathname1 '\' filename_sorted(i).name];
    load(fileDirName);

    idx_find =[];
    Light = [];
    idx_find = [];
    idx_pulse = [];
    Light = data(:,5);
    StimThr = 2; % set a threshold to filter light on period %TTR output 5 volt(actually ~4.5)
    idx_find = find(Light>StimThr);
    idx_pulse=find(diff(idx_find)>1);
   
   

    for pp = 1:5 % n to check
        ntrial =3;
        c_col = pp*2; % current column value
        if extract_guide(c_col) == -1
            continue       
        elseif pp == 1 
            ntrial = extract_guide(1);
            p_idx = 1;
            varlist = 4;
        elseif pp == 2 
            p_idx = 1;
            varlist = 2;
        elseif pp == 3 
            p_idx = 1;
            varlist = 1;
        elseif pp == 4 
            p_idx = 2;
            varlist = 3; 
        elseif pp == 5 
            p_idx = 3;
            varlist = 2;
        end
        extract_s = extract_guide(c_col);
        extract_e = extract_guide(c_col+1);

        %jrange_low = [];
        %jrange_high = [];
        
        on = protocol_set(p_idx,5);
        if pp ==5
            jrange_low = on*1000/SampleRate;
        else
            jrange_low = on*1000/SampleRate/300;
        end
        jrange_high = on*1000/SampleRate*1.5+1000;
    
        current_time = 0;
        
        idx_interval = [];
        idx_interval_low = [];
        idx_interval_high = [];
        idx_interval = find(diff(idx_find)>=jrange_low);
        %idx_interval_high = find(diff(idx_find)<=jrange_high(j));
        %idx_interval = intersect(idx_interval_high,idx_interval_low);
        if length(idx_interval<ntrial) 
            idx_interval(end+1)= length(idx_find);
        end
        k=1;
        kk=1;

        if varlist == 3 %%% Special case for Blood Pressure analysis 
          if strcmp(info{order,2}, 'C7') | strcmp(info{order,2}, 'C14')
            if strcmp(info{order,2}, 'C7')
                  ntrial =2;  %C7PVHE only use the last two trials for BP due to unstable baseline
            elseif strcmp(info{order,2}, 'C14')
                  ntrial =3;  %C14 PVHE only use the last three trials for BP due to unstable baseline
            end
            k=length(idx_interval);
            current_time = idx_interval(k);
            while kk<=ntrial & k>=1 
                if idx_interval(1)<current_time
                    endidx = idx_interval(k);
                    OnEnd = idx_find(endidx);
                    %OnEnd = idx_find(endidx)+(1000/protocol_set(p_idx,2)-protocol_set(p_idx,3))/SampleRate;
                    
                    Start = int32(OnEnd-protocol_set(p_idx,5)/4*1000/SampleRate*8+1);
                    %End = int32(OnEnd+protocol_set(p_idx,5)/4*1000/SampleRate);
                    
                    if pp==5 
                        End = int32(OnEnd+protocol_set(p_idx,5)*1000/SampleRate);
                    else
                        End = int32(OnEnd+protocol_set(p_idx,5)/20*1000/SampleRate);
                    end
                    OnStart = int32(OnEnd-protocol_set(p_idx,5)*1000/SampleRate+1);
                    
    
                    if OnStart<(extract_s*1000/SampleRate) | OnEnd>(extract_e*1000/SampleRate)
                        k=k-1;
                        display 'No0'
                        continue;
                    end
                   if Start<0 | OnStart<0
                         k = k-1;
                            display 'No1'
                        continue;
                    else
                        if sum(Light(Start:(OnStart-4000/protocol_set(p_idx,2)/SampleRate))>StimThr)>protocol_set(p_idx,3)
                            k = k-1;
                            display 'No2'
                            continue;
                        else
                            if sum(Light(OnStart:(OnStart+1000)))==0
                            %if sum(Light(OnStart:(OnStart+4000/protocol_set(p_idx,2)/SampleRate)))==0
                                k = k-1;
                                display 'No3'
                                continue;
                            end
                        end
                    end
                    display 'Yes'
                    
                    current_time = endidx;

                    trial_order(p_idx,varlist)= trial_order(p_idx,varlist) +1;
                    trialinfo_extract{p_idx,varlist}(trial_order(p_idx,varlist),:) = [info(order,:) kk];
                    data_extract{p_idx,varlist}(:, trial_order(p_idx,varlist))=data(Start:End,varlist);
                    
                    if pp==5 | pp==4 % Extract heart rate data with blood pressure
                        varlist_extra = 1;
                        trial_order(p_idx,varlist_extra)= trial_order(p_idx,varlist_extra) +1;
                        trialinfo_extract{p_idx,varlist_extra}(trial_order(p_idx,varlist_extra),:) = [info(order,:) kk];
                        data_extract{p_idx,varlist_extra}(:, trial_order(p_idx,varlist_extra))=data(Start:End,varlist_extra);
                    end
                    
                    kk = kk+1;
                    k = k-1;
                    
                else
                    k = k-1;
                end
            end
            continue;
          end
        end    
        current_time = 0;

        %%% Check if second loop is running
      %   disp(['Running second loop? varlist=', num2str(varlist), ...
      % ', info{order,2}= ', info{order,2}]);
      % 
      %   disp(['kk=', num2str(kk), ', ntrial=', num2str(ntrial), ...
      %         ', k=', num2str(k), ', len(idx_interval)=', num2str(length(idx_interval))]);
        while kk<=ntrial & k<=length(idx_interval) 
            if idx_interval(k)>current_time
                endidx = idx_interval(k);
                OnEnd = idx_find(endidx);
                %OnEnd = idx_find(endidx)+(1000/protocol_set(p_idx,2)-protocol_set(p_idx,3))/SampleRate;
                
                Start = int32(OnEnd-protocol_set(p_idx,5)/4*1000/SampleRate*8+1);
                %End = int32(OnEnd+protocol_set(p_idx,5)/4*1000/SampleRate);
                
                if pp==5 
                    End = int32(OnEnd+protocol_set(p_idx,5)*1000/SampleRate);
                else
                    End = int32(OnEnd+protocol_set(p_idx,5)/20*1000/SampleRate);
                end
                OnStart = int32(OnEnd-protocol_set(p_idx,5)*1000/SampleRate+1);
                

                if OnStart<(extract_s*1000/SampleRate) | OnEnd>(extract_e*1000/SampleRate)
                    k=k+1;
                    display 'No0'
                    continue
                end
               if Start<0 | OnStart<0
                     k = k+1;
                        display 'No1'
                    continue
                else
                    if sum(Light(Start:(OnStart-4000/protocol_set(p_idx,2)/SampleRate))>StimThr)>protocol_set(p_idx,3)
                        k = k+1;
                        display 'No2'
                        continue
                    else
                        if sum(Light(OnStart:(OnStart+1000)))==0
                        %if sum(Light(OnStart:(OnStart+4000/protocol_set(p_idx,2)/SampleRate)))==0
                            k = k+1;
                            display 'No3'
                            continue
                        end
                    end
                end
                display 'Yes'
                
                current_time = endidx;
                trial_order(p_idx,varlist)= trial_order(p_idx,varlist) +1;
                trialinfo_extract{p_idx,varlist}(trial_order(p_idx,varlist),:) = [info(order,:) kk];
                data_extract{p_idx,varlist}(:, trial_order(p_idx,varlist))=data(Start:End,varlist);
                if pp==5 | pp==4 % Extract heart rate data with breathing
                    varlist_extra = 1;
                    trial_order(p_idx,varlist_extra)= trial_order(p_idx,varlist_extra) +1;
                    trialinfo_extract{p_idx,varlist_extra}(trial_order(p_idx,varlist_extra),:) = [info(order,:) kk];
                    data_extract{p_idx,varlist_extra}(:, trial_order(p_idx,varlist_extra))=data(Start:End,varlist_extra);
                end
                kk = kk+1;
                k = k+1;
                
            else
                k=k+1;
            end
        end

 
    end

end

info_tb=cell2table(info,'VariableNames',{'id','name','group','type'});
