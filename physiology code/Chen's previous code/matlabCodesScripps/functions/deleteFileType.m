function temp = deleteFileType(fileDirectory, fileType)
%% delete all files of a certain kind inside a folder
% file typeexamples: 'tif','mat',etc. Not . in string
% Ran 20200101
    cd(fileDirectory);
    files = dir('*.' fileType);
    for iFile = 1:length(files)
        eval(['delete ' files(iFile).name]);
    end
end



