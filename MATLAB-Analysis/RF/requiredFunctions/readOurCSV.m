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