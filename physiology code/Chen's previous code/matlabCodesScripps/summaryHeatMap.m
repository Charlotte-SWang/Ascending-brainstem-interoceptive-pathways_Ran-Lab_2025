%% draw heat maps after pooling all responses together
% Ran 20190708
%% master file analysis
% Ran 20190127
% Pay attention to Lines 119-121 (not all ora stimuli have three trials)

clear all; close all; clc;
addpath('F:\matlabCodes\functions');
summaryDataDir = ['E:\tempAnalysis\202011 summaries\AllData'];

% blockList = {'ora','lar','lst','int','ili','cec'}; % (not all ora stimuli have three trials). Will show error in "allRes(cellCount,:) = dFFadj(:,iCell)"
% blockList = {'lar','lst','int'};
% blockList = {'PBS','Hla','CAc'};
% blockList = {'glc','nac','sal'}; % pay attention to Line 161 when analyzing int chem
% blockList = {'lar','lst','int'};
% blockList = {'ill','ili'};
% blockList = {'glc','int','lst'};
% blockList = {'glc','nac','sal'};
% blockList = {'glc','sal'};
% blockList = {'lst','lar'};
blockList = {'lst','lst'};
plotLength = 35; % plot 30 frames of response
plotResponder = 0;  % if 1, plot all responsive neurons. if 0, plot all non-responders (also Lines163-164)
ifSorting = 0; 
preStimFrame = 5; % number of frames before a stimulus to be included in the plot
redefineBlocks = 1;

fileName = dir([summaryDataDir '\*.mat']);
fileNumber = length(fileName);
cellCount = 0;
count = 0;
count2 = 0;
narrowTuneThd = 0;
cellNumberThd = 5;                     % if the number of total responders is smaller than cellNumberThd, skip the FOV
countSkipFile = 0;
countFile = 0;
nBlock = length(blockList);
for iFile = 1:fileNumber
    load([summaryDataDir '\' fileName(iFile).name],'blockName','datadir','maxdFF','goodCellBlock','datadirt','dFF','stimBlocks','isBadCell','F0Each','nCell','stim','maxdFFExt','thd','maxEachBlock','maxEachStimMinusThd','maxEachBMinusThd','maxEachBSelMinusThd','Fstd','goodCellEach','stimName','blockSto','nCell');
    goodEachIndiv{iFile} = goodCellEach;
    skipFile(iFile) = 0;
    clear block;

    %% use the following lines if only blocks are specified
    blockName(find(contains(blockName,'sto'))) = {'lst'};           % historical reason...
    blockName(find(contains(blockName,'HLa'))) = {'lar'};           % historical reason...
%     if str2num(datadirt(12:19)) < 20191126 | str2num(datadirt(12:19)) >= 20200210&str2num(datadirt(12:19)) < 20200822
%         blockName(find(contains(blockName,'ili'))) = {'wil'};       % find the ilium stim was too anterior and exclude them from the analysis
%     end
    timeWindow = [];
    indivWindow = [];
    count = 0;
    
    for iBlock = 1:nBlock
        oriBlockID = find(contains(blockName,blockList{iBlock}),1,'First');
        if isempty(oriBlockID) == 1 | length(stimBlocks{oriBlockID}) < 4 | oriBlockID ~= 1 | length(stim{1}) == 30
            skipFile(iFile) = 1;
            break;
        else
            block(iBlock) = oriBlockID;
        end        
    end
    if skipFile(iFile) == 1
        countSkipFile = countSkipFile + 1;
        continue;
    end  
                        %     sumGoodBlock = sum(goodCellBlock(:,block),1); % calculate the number of cells responding to each block so that experiments with fewer responders than the threshold would be excluded
                        % %         if any(sumGoodBlock(1:3)<cellNumberThd)    % int chem
                        % %     if any(sumGoodBlock(1:2)<cellNumberThd)    % int chem
                        %         if any(sumGoodBlock(1)<cellNumberThd)    % int chem
                        % %     if sumGoodBlock(end)<cellNumberThd    % use this line only when analyzing 50 uL int and ili data
                        %         skipFile(iFile) = 1;
                        %     end
                            
                            % nStim = count;
                            % if skipFile(iFile) == 1
                            %     countSkipFile = countSkipFile + 1;
                            %     continue;
                            % end
                            % temp = maxEachBMinusThd(:,block);
                            % veryMax = max(temp,[],2);
                            % relaResp = temp ./ veryMax;
                        %     
    regStim = [1 2 3 4];
    nStim = length(regStim);
    dFFadj=[];              % similar to the dFFadj in the other files
    dFFadjIndiv=[];
    for iStim = 1:nStim
        % renormalize the response using the F0 for each stim, using the first 5 frames after the stimulation onset or the background for each stimulus, whichever is higher
        junk = mean(dFF(stim{regStim(iStim)}(1:5),:),1);
        biggerF0 = max(junk,F0Each(regStim(iStim),:));
        biggerF0 = max(0,biggerF0);
        tempStim = stim{regStim(iStim)}(1) - preStimFrame : stim{regStim(iStim)}(1) + plotLength - 1 ;   
        if iStim == 1
             dFFadj = [dFFadj;dFF(tempStim,:)-biggerF0];
        else
            dFFadj = [dFFadj;dFF(tempStim,:)];
        end
        tempStimIndiv = stim{regStim(iStim)}(1) - preStimFrame : stim{regStim(iStim)}(1) + plotLength - 1; % include five more frames before stimulation onset
        dFFadjIndiv = [dFFadjIndiv;dFF(tempStimIndiv,:)-biggerF0];  
%         dFFadj = [dFFadj;dFF(stim{regStim(iStim)},:)-biggerF0];
    end
            
    for iCell = 1:nCell
        if any(goodCellEach(iCell,regStim),2)== 1 & isBadCell(iCell) == 0 
% %         if sum(goodCellEach(iCell,regStim),2)== 0 & isBadCell(iCell) == 0
%         if isBadCell(iCell) == 0
            if goodCellEach(iCell,2) == 0 & maxdFF(iCell,regStim(2))/maxdFF(iCell,regStim(3)) > 0.3
                continue;
            end
            cellCount = cellCount + 1;
%             alldFF(cellCount,:) = dFF(:,iCell);
            allRes(cellCount,:) = dFFadj(:,iCell);
            goodBlockAll(cellCount,:) = goodCellBlock(iCell, block);
            goodCellAll(cellCount,:) = goodCellEach(iCell,regStim); 
            maxEachBlockAll(cellCount,:) = maxEachBlock(iCell,block);
            maxdFFAll(cellCount,:) = maxdFF(iCell,regStim);
            maxEachBlockAllMinusThd(cellCount,:) = maxEachStimMinusThd(iCell,block);
            normRes(cellCount,:) = dFFadj(:,iCell)./max(maxdFFAll(cellCount,:));
            allIndivRes(cellCount,:) = dFFadjIndiv(:,iCell);
            if plotResponder == 1
                firstStimToResp = find(goodCellAll(cellCount,:) == 1,1,'First');  % find the first stimulus that evokes a positive response, and use the max amplitude to that stimulus for sorting
                maxdFFTo1stStim(cellCount) = maxdFF(iCell,regStim(firstStimToResp));
            end
        end
    end
    maxEachBlockAll(maxEachBlockAll<0)=0;
    count2 = count2 + 1;
    percentageResponder(count2,:) = sum(goodCellBlock(:, block),1)/nCell;
    nCellAll(count2) = nCell;
    allGoodCell(iFile) = sum(any(goodCellEach,2)==1);
    nBadCell(iFile) = sum(isBadCell);
    percentBadCell(iFile) = sum(isBadCell)/nCell;
    % count the number of animals included in this analysis
    idx = find(datadirt=='\',2,'first');junk = idx(end);
    allAnimalID(iFile,:) = datadirt(junk+1:junk+8);
    miceCount = size(unique(allAnimalID,'rows'),1); 
end
% 
%% calculate correlation matrix

resList{1} = find(goodCellAll(:,2) == 1);
resList{2} = find(goodCellAll(:,2) == 0);


% resList{1} = find(goodBlockAll(:,1) == 1 & goodBlockAll(:,2) == 1 & goodBlockAll(:,4) == 1);
% resList{2} = find(goodBlockAll(:,1) == 1 & goodBlockAll(:,2) == 1 & goodBlockAll(:,4) == 0);
% resList{3} = find(goodBlockAll(:,1) == 1 & goodBlockAll(:,2) == 0 & goodBlockAll(:,4) == 1);
% resList{4} = find(goodBlockAll(:,1) == 0 & goodBlockAll(:,2) == 1 & goodBlockAll(:,4) == 1);
% resList{5} = find(goodBlockAll(:,1) == 1 & goodBlockAll(:,2) == 0 & goodBlockAll(:,4) == 0);
% resList{6} = find(goodBlockAll(:,1) == 0 & goodBlockAll(:,2) == 1 & goodBlockAll(:,4) == 0);
% resList{7} = find(goodBlockAll(:,1) == 0 & goodBlockAll(:,2) == 0 & goodBlockAll(:,4) == 1);


% resList{1} = find(any(goodBlockAll(:,1:3),2) == 1 & goodBlockAll(:,4) == 1);
% resList{2} = find(any(goodBlockAll(:,1:3),2) == 1 & goodBlockAll(:,4) == 0);
% resList{3} = find(any(goodBlockAll(:,1:3),2) == 0 & goodBlockAll(:,4) == 1);
% resList{4} = find(any(goodBlockAll(:,1:3),2) == 0 & goodBlockAll(:,4) == 0 & goodBlockAll(:,5) == 1);
% % resList{1} = find(any(goodBlockAll,2) == 1);
% resList{1} = find(sum(goodBlockAll,2) == 0);
% 

% % lower vs upper stom
% aa = maxdFFAll(:,1) ./ maxdFFAll(:,3);
% resList{1} = find(aa >= 0.7);
% resList{2} = find(aa < 0.7);
% % linear regression
% volume = [0.15;0.3;0.6;0.9];
% for iCell = 1:cellCount
%     junk = maxdFFAll(iCell,5:8)';
%     b1(iCell) = volume\junk;
% end
% resList{1} = find(b1 > median(b1));
% resList{2} = find(b1 < median(b1));
%     


% % create a helper matrix a to go through all combinations of 6 stim
% a = zeros(2,nBlock); a(2,:) = 1; 
% allcombmatrix = allcomb(a(:,1), a(:,2), a(:,3), a(:,4), a(:,5), a(:,6)); % this line needs to be modified if more than 6 stim  , a(:,4), a(:,5), a(:,6)
% allcombmatrix = allcombmatrix(2:end,:); % use this line to exclude cells that do not respond to any stim
% allcombmatrix = flip(allcombmatrix,2);
% % use the following three lines to sort the allcombmatrix so that broadly-tuned cells are at the bottom of the matrix
% temp = sum(allcombmatrix,2);
% [junk,I] = sort(temp);
% allcombmatrix = allcombmatrix(I,:);
% for iRespProfile = 1:size(allcombmatrix,1)
%     resList{iRespProfile} = find(sum(goodBlockAll == allcombmatrix(iRespProfile,:),2) == nBlock);
% %     [~,resList{iRespProfile}] = intersect(goodBlockAll,allcombmatrix(iRespProfile,:),'rows');
%     lenResList(iRespProfile) = length(resList{iRespProfile});
% end
% 


plotList = [];
temp1 = mean(maxdFFAll,2);
% for iList = 1:length(resList)
for iList = 1:length(resList)
    lenList(iList) = length(resList{iList});
    if ifSorting == 1
%         [junk,temp]=sort(maxdFFTo1stStim(resList{iList}),'descend');
        [junk,temp]=sort(temp1(resList{iList}),'descend');
        plotList = [plotList;resList{iList}(temp)];
    else
        plotList = [plotList;resList{iList}];
    end
end
% 
% % sorting according to the relative response to stom vs int to put most
% % selective cells at the top and bottom
% % plotList = 1:cellCount;
% stoVInt = maxEachBlockAll(:,1) ./ maxEachBlockAll(:,2);
% [junk2,temp2]=sort(stoVInt,'descend');
% figure;
% hxx1 = imagesc(normRes(temp2,:),[0.3 1]);
% % hxx1 = imagesc(normRes(temp2(7820:7919),:),[0.3 1]);
% colormap hot;
% 
% % % use the following lines to average responses across repeats
% % temp3(:,:,1) = allRes(:,[1:50 101:150 201:250]);
% % temp3(:,:,2) = allRes(:,[51:100 151:200 251:300]);
% % temp2 = mean(temp3,3);

figure;
% hxx1 = imagesc(normRes(plotList,:),[0.1 1]); colormap hot;
hxx1 = imagesc(allRes(plotList,:),[10 70]); colormap jet;
% hxx1 = imagesc(temp2(plotList,:),[10 100]); colormap jet;


for iList = 1:length(resList)
    for iBlock = 1:length(blockList)
        a(iBlock,(iList-1)*3+1) = mean(maxEachBlockAll(resList{iList},iBlock));
        a(iBlock,(iList-1)*3+2) = std(maxEachBlockAll(resList{iList},iBlock));
        a(iBlock,(iList-1)*3+3) = lenList(iList);
    end
end
% hxx1 = imagesc(normRes(plotList,:),[0.1 1]);
% colormap hot;
% title(sprintf('Res to sto, int, NaCl'));
% saveas(gcf,[summaryDataDir 'Stomach, intestinal distension and chemical responses.jpg']);
% saveas(gcf,[summaryDataDir 'Stomach, intestinal distension and chemical responses.eps'],'psc2');


