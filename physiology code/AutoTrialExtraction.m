%% Extract trial information
% Shiqi Wang 20250114

clear all


%% Input

%%% Attention !!!%%%
SampleRate = 0.2; % default unit: ms --check 'isi' & 'isi_units'

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
filename_sorted=dir('*_*_*_*_*.mat');
n=length(filename_sorted);%file numbers
data_extract = {};
for i = 1:n
    % Prepare information
    [~,token,~] = fileparts(filename_sorted(i).name);
    token_split = split(token,'_');

    order = str2num(token_split{1});


    info{order,1} = order;
    info{order,2} = token_split{2};
    info{order,3}= token_split{3};
    time_split = split(token_split{4},'-');
    info{order,4} = string(join({time_split{1},time_split{2}},'-'));
    info{order,5} = time_split{3};

    % Extract preparation
    extract_guide = [];
    [~, rowNum] = ismember(filename_sorted(i).name, individual_protocol.Properties.RowNames);
    extract_guide = individual_protocol{rowNum,:};
    m = max(extract_guide);


    % Read data file

    % Read data
    fileDirName = [pathname1 '\' filename_sorted(i).name];
    load(fileDirName);
  
    idx_find =[];
    Light = [];
    idx_find = [];
    idx_pulse = [];
    Light = data(:,5);
    StimThr = 1; % set a threshold to filter light on period %TTR output 5 volt(actually ~4.5)
    idx_find = find(Light>StimThr);
    idx_pulse=find(diff(idx_find)>1);
   

    jrange_low = [];
    jrange_high = [];
    for j =1:m
        p_idx = find(extract_guide==j);
        on = protocol_set(p_idx,4);
        jrange_low(j) = on*1000/SampleRate-10;
        jrange_high(j) = on*1000/SampleRate*1.5+10;
    end
    
    current_time = 0;

    for j = 1:m
        p_idx = find(extract_guide==j);
        idx_interval = [];
        idx_interval_low = [];
        idx_interval_high = [];
        idx_interval = find(diff(idx_find)>=jrange_low(j));
        %idx_interval_high = find(diff(idx_find)<=jrange_high(j));
        %idx_interval = intersect(idx_interval_high,idx_interval_low);
        if length(idx_interval<3)
            idx_interval(end+1)= length(idx_find);
        end
        k=1;
        kk=1;

        while kk<=3 & k<=length(idx_interval)
            if idx_interval(k)>current_time
                endidx = idx_interval(k);
                OnEnd = idx_find(endidx);
                
                Start = int32(OnEnd-protocol_set(p_idx,4)/4*1000/SampleRate*8+1);
                End = int32(OnEnd+protocol_set(p_idx,4)/4*1000/SampleRate);
                
                OnStart = int32(OnEnd-protocol_set(p_idx,4)*1000/SampleRate+1);
                if Start<0 | OnStart<0
                     k = k+1;
                        display 'No'
                    continue
                else
                    if sum(Light(Start:(OnStart-1000/protocol_set(p_idx,2)/SampleRate))>StimThr)>protocol_set(p_idx,3)
                        k = k+1;
                        display 'No'
                        continue
                    else
                        if sum(Light(OnStart:(OnStart+1000/protocol_set(p_idx,2)/SampleRate)))==0
                            k = k+1;
                            display 'No'
                            continue
                        end
                    end
                end
                display 'Yes'
                
                current_time = endidx;
                var_list = find(protocol_set(p_idx,(FixInfo+1):end)==1);
                for l = 1:length(var_list)
                    trial_order=3*(order-1)+kk;
                    data_extract{p_idx,var_list(l)}(:,trial_order)=data(Start:End,var_list(l));
                end
                         
                kk = kk+1;
                k = k+1;
                
            else
                k=k+1;
            end
        end

    end

end
info_tb=cell2table(info,'VariableNames',{'id','group','name','date','time(hours)'});
