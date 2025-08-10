function [block,timeWindow,listStim] = timeWindowInBlocks(stimBlocks,blockName,listBlock,stim,plotLengthEach,preStimFrame)
% Determine the time window in blocks specified in listBlock
% Ran 20201104


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

end

