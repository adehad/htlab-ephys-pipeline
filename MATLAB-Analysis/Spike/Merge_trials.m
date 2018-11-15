%%%% Merge_trials
%%%% A. Haddad [Sep 2018]

%% Merges .bin files
% Specify the file names in namePart1 and numFiles
% also saves mergeOutput.txt, that list the length (in samples) of each
% individual file
% Set paramters of data
nChans = 4;
fSamp = 30e3; % not actually needed for script

% fileNameList = [ '131210_7986_MovingObjects_1.bin' ...
%                 ];

numFiles = [1]; % note: an array
namePart1 = 'Wing_ephys_181028_converted';
namePart2 = '.bin';
for k=1:length(numFiles)
    fileNameList(k,:) = [ namePart1, num2str(numFiles(k)), namePart2 ];
end
            
fileSizeList = zeros(1,size(fileNameList,1)); % patch zeros

% Open file(s) & concatenate
dataRAW = [];
for kk=1:size(fileNameList,1)            
    fileID = fopen(fileNameList(kk,:));

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

% readOurCSV function (defined at the bottom of the script)
% to read the above created file and store useful info into struct
outputStruct = readOurCSV(fileName)



% save merged
fileID = fopen([namePart1(1:6),'_merged46.bin'],'w');
fwrite(fileID, dataRAW, 'int16','l'); % little endian write
fclose(fileID);

clear temp

function outputStruct = readOurCSV(theFileName)
% Parses the merge_info text file to bring useful information into the
% matlab workspace
%   Input a filename that can be opened by fopen - e.g. char vector
%   The text is comma separated of the form:
%   Text [File name.bin], is, Samples per .bin, samples long
%   The 'filename' and 'samples' are stored in the output struct
%   additionally an 'index' is created that stores the location in the
%   merged bin file

fileID = fopen(theFileName,'r');
% Format for each line of text:
%   column1: text (%s)
%	column2: categorical (%C)
%   column3: double (%f)
%	column4: categorical (%C)
formatSpec = '%s%C%f%C%[^\n\r]';
text = textscan(fileID,formatSpec,'Delimiter',',');
fclose(fileID);

% Store filenames
for inc = 1:size(text{1,1},1)
    outputStruct.filename(inc,:) = char(text{1,1}(inc));
end

% Store Sample Sizes
for inc = 1:size(text{1,3},1)
    outputStruct.samples(inc,:) = double(text{1,3}(inc));
end

% Store Sample Index (relative to location in the merged .bin file)
tot = 0;
for inc = 1:size(text{1,3},1)
    tot = outputStruct.samples(inc) + tot;
    outputStruct.index(inc,:) = tot;
end

end