%% Extract trial information
% Shiqi Wang 20250114

%% Input
%Read instructions for protocol extraction
[FileName,PathName]=uigetfile('*.csv');
addpath(PathName)
protocol_preview = readtable(FileName)
protocol = readmatrix(FileName)
FixInfo = 6 % Channel is from 7:end


PulseL = 0.01; % pulse length, unit:s
PulseF = 20; % pulse frequency, unit:Hz
PulseD = 10; % Duration of each group of pulse
   
%%% Attention !!!%%%
SampleRate = 0.2; % default unit: ms --check 'isi' & 'isi_units'
   
   for i = 1:nStim
        OnStart(i) = idx_find(ii);
        OnDuration(i) = Stim(i)*1000/SampleRate;
        BaseLength(i) = BaseStim(i)*1000/SampleRate;
        OnEnd(i) = OnStart(i)+OnDuration(i)-1; 
        BaselineStart(i) = OnStart(i)-BaseLength(i);
        BaselineEnd(i) = OnStart(i)-1;
        OffStart(i) = OnEnd(i)+1;
        ii = find(idx_find>OffStart(i),1);
 
        ECG_analysis{i} =ECG(BaselineStart(i):OnEnd(i)) ;
        BP_analysis{i} =BP(BaselineStart(i):OnEnd(i)) ;
        Breath_analysis{i} =Breath(BaselineStart(i):OnEnd(i)) ;
        Light_analysis{i}= Light(BaselineStart(i):OnEnd(i));
        idx{i} = find(Light_analysis{i}>StimThr);
        id{i} = '10s on off'
    end  
%% 


timeWindow = [];            % the time points that belong to the blocks
nStim = 0;
listStim = [];

for iBlock = 1:length(listBlock)
    oriBlockID = find(contains(blockName,listBlock{iBlock}));
    if isempty(oriBlockID) ~= 1
        block(iBlock) = oriBlockID;
        for iStim = 1:length(stimBlocks{oriBlockID})
            if nargin < 6
                timeWindow = [timeWindow stim{stimBlocks{oriBlockID}(iStim)}(1:plotLengthEach)];
            end
        end
        listStim = [listStim stimBlocks{oriBlockID}];
    else
        block(iBlock) = NaN;
    end
end

fil