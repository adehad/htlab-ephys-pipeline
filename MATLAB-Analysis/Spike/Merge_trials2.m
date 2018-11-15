%%%% Merge_trials
%%%% A. Haddad [Sep 2018]

%% Merges .bin files
% Specify the file names in namePart1 and numFiles
% also saves mergeOutput.txt, that list the length (in samples) of each
% individual file
% Set paramters of data
nChans = 2;
fSamp = 40000; % not actually needed for script

%% conpile file names
namePart1 = '131217_8014p_MovingObjects_'; %'131222_8021_MovingObjects_'; %'131210_7986_MovingObjects_';131228_8078_MovingObjects_  '131230_DE10_MovingObjects_'
numFiles = [1:7 9]; % note: an array
namePart2 = '.bin';
for k=1:length(numFiles)
    fileNameList{k,:} = [ namePart1, num2str(numFiles(k)), namePart2 ];
end
    %repeat again to add more files
lSize=size(fileNameList,1); % first list size
namePart1 = '131217_8014q_MovingObjects_'; %'131222_8021q_MovingObjects_'; %'131210_7986_MovingObjects_';131228_8078_MovingObjects_  '131230_DE10_MovingObjects_'
numFiles = [1:2]; % note: an array
namePart2 = '.bin';
for k=1:length(numFiles)
    fileNameList{k+lSize,:} = [ namePart1, num2str(numFiles(k)), namePart2 ];
end

fileSizeList = zeros(1,size(fileNameList,1)); % patch zeros

%% Open file(s) & concatenate
dataRAW = [];
for kk=1:size(fileNameList,1)            
    fileID = fopen(fileNameList{kk,:});

    temp = fread(fileID,[nChans, Inf], 'int16','l'); % little endian open
    dataRAW = cat(2,dataRAW,temp);

    fileSizeList(kk) = length(temp);
    fclose(fileID);
end
dataRAW = int16(dataRAW);

% save details of merge

% % CSV file
% fileName = [namePart1(1:6),'_(',char(strjoin(string(numFiles),'+')),')_merge_info.csv'];
fileName = [namePart1(1:6),'_merge_info.csv'];
fileID = fopen(fileName,'w');
formatSpec = '%s, is ,%.0f, samples long\n';
fprintf(fileID,formatSpec,[string(fileNameList)';string(fileSizeList)]);
fclose(fileID);

%% save merged
fileID = fopen([namePart1(1:6),'_merged.bin'],'w');
fwrite(fileID, dataRAW, 'int16','l'); % little endian write
fclose(fileID);
clear temp