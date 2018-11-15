%%%% Merge_trials
%%%% A. Haddad & HT. Lin [Sep 2018]
%%%% This script merges .bin files
% Specify the file names in namePart1 and numFiles
% also saves mergeOutput.txt, that list the length (in samples) of each
% individual file
% Set paramters of data

%% define recording paramters
nChans = 5; %2
fSamp = 40000; % not actually needed for script
% for data collected via SpikeGL, we can incorporate MetafileRead to
% extract this information

%% conpile file names
namePart1 = '150526__MovingObjects_'; %'131210_7986q_MovingObjects_'; '131217_8014p_MovingObjects_'; %'131222_8021_MovingObjects_'; %'131210_7986_MovingObjects_';131228_8078_MovingObjects_  '131230_DE10_MovingObjects_'
numFiles = [1:10]; % note: an array
namePart2 = '.bin';
for k=1:length(numFiles)
    fileNameList{k,:} = [ namePart1, num2str(numFiles(k)), namePart2 ];
end

%     %repeat again to add more files
% lSize=size(fileNameList,1); % first list size
% namePart1 = '131217_8014q_MovingObjects_'; %'131222_8021q_MovingObjects_'; %'131210_7986_MovingObjects_';131228_8078_MovingObjects_  '131230_DE10_MovingObjects_'
% numFiles = [1:2]; % note: an array
% namePart2 = '.bin';
% for k=1:length(numFiles)
%     fileNameList{k+lSize,:} = [ namePart1, num2str(numFiles(k)), namePart2 ];
% end

fileSizeList = zeros(1,size(fileNameList,1)); % patch zeros

%% open file(s) & concatenate
dataRAW = [];
for kk=1:size(fileNameList,1)            
    fileID = fopen(fileNameList{kk,:});

    temp = fread(fileID,[nChans, Inf], 'int16','l'); % little endian open
    fileSizeList(kk) = length(temp);
    fclose(fileID);      
%     dataRAW = cat(2,dataRAW,temp);

    % save merged
%     fileID = fopen([namePart1(1:6),'_merged.bin'],'w');
    fileID = fopen([namePart1(1:6),'_merged.bin'],'a');
    fwrite(fileID, temp, 'int16','l'); % little endian write
    fclose(fileID);
    clear temp
end
dataRAW = int16(dataRAW);

%% save details of the merge operation
% CSV file
% fileName = [namePart1(1:6),'_(',char(strjoin(string(numFiles),'+')),')_merge_info.csv'];
fileName = [namePart1(1:6),'_merge_info.csv'];
fileID = fopen(fileName,'w');
formatSpec = '%s, is ,%.0f, samples long\n';
fprintf(fileID,formatSpec,[string(fileNameList)';string(fileSizeList)]);
fclose(fileID);
