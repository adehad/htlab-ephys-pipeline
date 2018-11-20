function binToBin(pathToDataFolder,dataNames,overwriteFiles, ...
    dataCh,nChDesired,interlaceCh,invertCh,mergeCh)
%% Saves OPEN-EPHYS Data as Binary & Interlaces, Adds Dummy channels, Saves ADC
% Last Updated: 11/11/2018
%
% INPUT:
%       pathToDataFolder (path to OpenEphys Data folders (dataFolderNames).
%           Uses current directory if argument is empty)
%       dataNames (string array of data .bin to
%           merge, 'all' chooses all folders in pathToDataFolder or current
%           dir)
%       dataCh          (Number of data channels in input .bin)
%       nChDesired      (if greater than length(dataChan) will pad with
%                          dummy channels of 1s)
%       interlaceCh     (if you want to interlace channel data)
%       invertCh        (1 if you want to invert)
%       mergeCh         (1 if you want merged data output)
%
% OUTPUT:                 Stored in the same folder as dataFolderNames
%      dataFolderName.bin (will groups channels to a single binary file)
%                           each row is a different channel
%      dataFolderName_converted.bin (above, but with dummy channels)
%                                    each row is a different channel
%      dataFolderName_ADC.bin        (ADC channels, never padded)
%      dataFolderName_timeStamps.mat (timestamp information - extracted from ADC )
%
% Functions required at the end of the file: none

% shorthand functions
isAnInteger = @(x) isfinite(x) & x == floor(x);
char2str = @(x) convertCharsToStrings(x);

% -----------------------------------------------
%                   INPUT
% -----------------------------------------------
% Find folders if using all folders in path
if strcmp(dataNames, 'all') || strcmp(dataNames, 'ALL')
    if ~isempty(pathToDataFolder)
        filesInPath = dir(pathToDataFolder);
    else
        filesInPath = dir;
    end
    folderFlags = [filesInPath.isdir] & ~strcmp({filesInPath.name},'.') & ~strcmp({filesInPath.name},'..');
    dataFolders = filesInPath(folderFlags);
    dataNames = char2str({dataFolders.name})';
end

% output file naming
namePart1 = erase(dataNames,filesep);
namePart1 = erase(namePart1,'.bin');
namePart2 = '.bin';

newName             = namePart1 + '_padded' + namePart2;
newNameINTERLACED   = namePart1 + '_interlaced' + namePart2;
newNameMERGED       = ['data_merged', namePart2];
newNameMETA         = 'merge_info.csv';

% if user is specifying a path to the data folder - i.e. not current folder
if ~isempty(pathToDataFolder)
    if pathToDataFolder(end) ~= filesep
        pathToDataFolder = [pathToDataFolder, filesep];
    end
    newName             = pathToDataFolder + newName;
    newNameINTERLACED   = pathToDataFolder + newNameINTERLACED;
    newNameMERGED       = [pathToDataFolder, newNameMERGED];
    newNameMETA         = [pathToDataFolder, newNameMETA];
end

% data channel extraction
if ~isempty(dataCh)
    % allows appending of merging data
    mergeLock = 1;
    
    for ii=1:length(dataNames)
        clear tempData
        
        % load data
        fileID = fopen(dataNames(ii));
        tempData = fread(fileID,[dataCh, Inf], 'int16','l'); % little endian open
        tempData = int16(tempData);
        fclose(fileID);
        
        % interlacing, if the option was enabled
        if interlaceCh
            if isAnInteger(dataCh/2) % is it an integer
                tempData = zeros(dataCh/2, size(tempData,2)*2 );
                for kk=1:dataCh/2
                    tempData(kk,1:2:end) = tempData(kk,:);
                    tempData(kk,2:2:end) = tempData(kk+(dataCh/2),:);
                end
            else
                warning(['InterlaceCh set to 1, but not enough channels to'...
                    'interlace. Skipping interlacing...'])
            end
            
            if overwriteFiles || ~isfile(newNameINTERLACED(ii))
                fileID = fopen(newNameINTERLACED(ii),'w');
                fileLock = 0;
            else
                fileLock = 1;
            end
        else
            if overwriteFiles || ~isfile(newName(ii))
                fileID = fopen(newName(ii),'w');
                fileLock = 0;
            else
                fileLock = 1;
            end
        end
        
        % create new split data .bin file
        formattedData = ones(nChDesired,length(tempData),'int16');
        formattedData(1:size(tempData,1),:) = tempData*(-1)^invertCh;
        formattedData = int16(formattedData);

        % save new split data .bin file, check for overwrites
        if ~fileLock
            fwrite(fileID, formattedData, 'int16','l'); % little endian write
            fclose(fileID);
            if fileID < 0
                warning('The split data could not be saved')
            end
        else
            warning('File exists and will not be overwritten');
        end
        
        % save merged data, check for overwrites
        if mergeCh
            % get size of each split
            fileSizeList(ii) = length(formattedData);
            
            if ~mergeLock
                fileID = fopen(newNameMERGED,'a');
            elseif isfile(newNameMERGED)
                if overwriteFiles
                    fileID = fopen(newNameMERGED,'w');
                    mergeLock = 0;
                else
                    warning('File exists and will not be overwritten');
                    mergeLock = 1;
                end
            else
                fileID = fopen(newNameMERGED,'w');
                mergeLock = 0;
            end
            
            if ~mergeLock
                fwrite(fileID, formattedData, 'int16','l');
                fclose(fileID);
                if fileID < 0
                    warning('The merged data could not be saved')
                end
            else
                warning('File exists and will not be overwritten');
            end
        end
    end
    
    % save merging metadata, check for overwrites
    if ~mergeLock
        fileID = fopen(newNameMETA,'w');
        formatSpec = '%s, is ,%.0f, samples long\n';
        fprintf(fileID,formatSpec,[namePart1'; string(fileSizeList)]);
        fclose(fileID);
    elseif mergeCh
        warning('File exists and will not be overwritten');
    end
else
    warning('Empty dataCh')
end

end