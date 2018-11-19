function OEtoBin(pathToDataFolder,dataFolderNames,overwriteFiles, ...
    dataCh,adcCh,nChDesired,interlaceCh,invertCh,mergeCh)
%% Saves OPEN-EPHYS Data as Binary & Interlaces, Adds Dummy channels, Saves ADC
% Last Updated: 11/11/2018
%
% INPUT:
%       pathToDataFolder (path to OpenEphys Data folders (dataFolderNames).
%           Uses current directory if argument is empty)
%       dataFolderNames (string array of folders with OpenEphys data to
%           merge, 'all' chooses all folders in pathToDataFolder or current
%           dir)
%       dataCh          (Number array that corresponds to : CH<> files
%                                                  leave empty if none)
%       adcCh           (Number array that corresponds to : ADC<> files
%                                                  leave empty if none)
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
% Functions required at the end of the file: load_open_ephys_data() & filesize()

% shorthand functions
isAnInteger = @(x) isfinite(x) & x == floor(x);
str2char = @(x) convertStringsToChars(x);
char2str = @(x) convertCharsToStrings(x);

% -----------------------------------------------
%                   INPUT
% -----------------------------------------------
% Find folders if using all folders in path
if strcmp(dataFolderNames, 'all') || strcmp(dataFolderNames, 'ALL')
    if ~isempty(pathToDataFolder)
        filesInPath = dir(pathToDataFolder);
    else
        filesInPath = dir;
    end
    folderFlags = [filesInPath.isdir] & ~strcmp({filesInPath.name},'.') & ~strcmp({filesInPath.name},'..');
    dataFolders = filesInPath(folderFlags);
    dataFolderNames = char2str({dataFolders.name})';
end

% output file naming
namePart1 = erase(dataFolderNames,filesep);
dataFolderNames = namePart1 + filesep;
namePart2 = '.bin';

newName             = dataFolderNames + namePart1 + '_padded' + namePart2;
newNameINTERLACED   = dataFolderNames + namePart1 + '_interlaced' + namePart2;
newNameADC          = dataFolderNames + namePart1 + '_ADC' + namePart2;
newNameTS           = dataFolderNames + namePart1 + '_timestamps.mat';
newNameMERGED       = ['data_merged', namePart2];
newNameMETA         = 'merge_info.csv';

% if user is specifying a path to the data folder - i.e. not current folder
if ~isempty(pathToDataFolder)
    if pathToDataFolder(end) ~= filesep
        pathToDataFolder = [pathToDataFolder, filesep];
    end
    newName             = pathToDataFolder + newName;
    newNameINTERLACED   = pathToDataFolder + newNameINTERLACED;
    newNameADC          = pathToDataFolder + newNameADC;
    newNameTS           = pathToDataFolder + newNameTS;
    newNameMERGED       = [pathToDataFolder, newNameMERGED];
    newNameMETA         = [pathToDataFolder, newNameMETA];
end

% data channel extraction
if ~isempty(dataCh)
    % allows appending of merging data
    mergeLock = 1;
    
    for ii=1:length(dataFolderNames)
        clear tempData
        for jj=1:length(dataCh)
            % load data
            [tempData(jj,:), ~, ~] = load_open_ephys_data([pathToDataFolder,...
                str2char(dataFolderNames(ii)),'100_CH',num2str(dataCh(jj)),'.continuous']);
            %disp(['finished loading Ch : ', num2str(dataCh(jj))]);
        end
        tempData = int16(tempData);
        
        % interlacing, if the option was enabled
        if interlaceCh
            if isAnInteger(length(dataCh)/2) % is it an integer
                tempData = zeros(length(dataCh)/2, size(tempData,2)*2 );
                for kk=1:length(dataCh)/2
                    tempData(kk,1:2:end) = tempData(kk,:);
                    tempData(kk,2:2:end) = tempData(kk+(length(dataCh)/2),:);
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

% adc channel extraction
if ~isempty(adcCh)
    for ii=1:length(dataFolderNames)
        clear tempADC
        for jj=1:length(adcCh)
            % load adc data
            [tempADC(jj,:), timeStamps, ~] = load_open_ephys_data([pathToDataFolder,...
                    str2char(dataFolderNames(ii)),'100_ADC',num2str(adcCh(jj)),'.continuous']);
            %disp(['finished loading ADC : ', num2str(adcCh(jj))]);
        end
        
        % mean(,2) column mean - in case of multiple channels where only one is bad
        if sum( mean( (tempADC) - double(int16(tempADC)) ,2) > 0.001 )
            warning(['The variation of the data is too small for integer ', ...
                'quantisation, mean quantisation error is: ', num2str(mean( (tempADC) - double(int16(tempADC)) )), ...
                newline, 'Multiplying data by 1000...']);
            tempADC = tempADC*1e3;
        end
        tempADC = int16(tempADC);

        % create new adc .bin file - Note: Does not get inverted
        formattedADC = zeros(length(adcCh), length(tempADC),'int16');
        formattedADC(1:size(tempADC,1),:) = tempADC(:,:);
        formattedADC = int16(formattedADC);

        % save new adc .bin file and timestamps, check for overwrites
        if overwriteFiles || ~isfile(newNameADC(ii))
            fileID = fopen(newNameADC(ii),'w');
            fileLock = 0;
        else
            fileLock = 1;
        end
        
        if ~fileLock
            save(newNameTS(ii),'timeStamps')
            fwrite(fileID, formattedADC, 'int16','l'); % little endian write
            fclose(fileID);
            if fileID > 0
                disp('Successfully saved new ADC data');
            else
                disp('The file could not be saved, check if the newFileNameADC is legal')
            end
        else
            warning('File exists and will not be overwritten');
        end
    end
else
    warning('Empty adcCh')
end
end

%%%% OPEN EPHYS DATA LOAD FUNCTION
%%%% https://github.com/open-ephys/analysis-tools/blob/master/load_open_ephys_data.m
function [data, timestamps, info] = load_open_ephys_data(filename)

%
% [data, timestamps, info] = load_open_ephys_data(filename)
%
%   Loads continuous, event, or spike data files into Matlab.
%
%   Inputs:
%
%     filename: path to file
%
%
%   Outputs:
%
%     data: either an array continuous samples (in microvolts),
%           a matrix of spike waveforms (in microvolts),
%           or an array of event channels (integers)
%
%     timestamps: in seconds
%
%     info: structure with header and other information
%
%
%
%   DISCLAIMER:
%
%   Both the Open Ephys data format and this m-file are works in progress.
%   There's no guarantee that they will preserve the integrity of your
%   data. They will both be updated rather frequently, so try to use the
%   most recent version of this file, if possible.
%
%

%
%     ------------------------------------------------------------------
%
%     Copyright (C) 2014 Open Ephys
%
%     ------------------------------------------------------------------
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     <http://www.gnu.org/licenses/>.
%

filetype = filename(max(strfind(filename,'.'))+1:end); % parse filetype

fid = fopen(filename);
filesize = getfilesize(fid);

% constants
NUM_HEADER_BYTES = 1024;
SAMPLES_PER_RECORD = 1024;
RECORD_MARKER = [0 1 2 3 4 5 6 7 8 255]';
RECORD_MARKER_V0 = [0 0 0 0 0 0 0 0 0 255]';

% constants for pre-allocating matrices:
MAX_NUMBER_OF_SPIKES = 1e6;
MAX_NUMBER_OF_RECORDS = 1e6;
MAX_NUMBER_OF_CONTINUOUS_SAMPLES = 1e8;
MAX_NUMBER_OF_EVENTS = 1e6;
SPIKE_PREALLOC_INTERVAL = 1e6;

%-----------------------------------------------------------------------
%------------------------- EVENT DATA ----------------------------------
%-----------------------------------------------------------------------

if strcmp(filetype, 'events')
    
    disp('Loading events file...');
    
    index = 0;
    
    hdr = fread(fid, NUM_HEADER_BYTES, 'char*1');
    eval(char(hdr'));
    info.header = header;
    
    if (isfield(info.header, 'version'))
        version = info.header.version;
    else
        version = 0.0;
    end
    
    % pre-allocate space for event data
    data = zeros(MAX_NUMBER_OF_EVENTS, 1);
    timestamps = zeros(MAX_NUMBER_OF_EVENTS, 1);
    info.sampleNum = zeros(MAX_NUMBER_OF_EVENTS, 1);
    info.nodeId = zeros(MAX_NUMBER_OF_EVENTS, 1);
    info.eventType = zeros(MAX_NUMBER_OF_EVENTS, 1);
    info.eventId = zeros(MAX_NUMBER_OF_EVENTS, 1);
    
    if (version >= 0.2)
        recordOffset = 15;
    else
        recordOffset = 13;
    end
    
    while ftell(fid) + recordOffset < filesize % at least one record remains
        
        index = index + 1;
        
        if (version >= 0.1)
            timestamps(index) = fread(fid, 1, 'int64', 0, 'l');
        else
            timestamps(index) = fread(fid, 1, 'uint64', 0, 'l');
        end
        
        
        info.sampleNum(index) = fread(fid, 1, 'int16'); % implemented after 11/16/12
        info.eventType(index) = fread(fid, 1, 'uint8');
        info.nodeId(index) = fread(fid, 1, 'uint8');
        info.eventId(index) = fread(fid, 1, 'uint8');
        data(index) = fread(fid, 1, 'uint8'); % save event channel as 'data' (maybe not the best thing to do)
        
        if version >= 0.2
            info.recordingNumber(index) = fread(fid, 1, 'uint16');
        end
        
    end
    
    % crop the arrays to the correct size
    data = data(1:index);
    timestamps = timestamps(1:index);
    info.sampleNum = info.sampleNum(1:index);
    info.nodeId = info.nodeId(1:index);
    info.eventType = info.eventType(1:index);
    info.eventId = info.eventId(1:index);
    
    %-----------------------------------------------------------------------
    %---------------------- CONTINUOUS DATA --------------------------------
    %-----------------------------------------------------------------------
    
elseif strcmp(filetype, 'continuous')
    
    disp(['Loading ' filename '...']);
    
    index = 0;
    
    hdr = fread(fid, NUM_HEADER_BYTES, 'char*1');
    eval(char(hdr'));
    info.header = header;
    
    if (isfield(info.header, 'version'))
        version = info.header.version;
    else
        version = 0.0;
    end
    
    % pre-allocate space for continuous data
    data = zeros(MAX_NUMBER_OF_CONTINUOUS_SAMPLES, 1);
    info.ts = zeros(1, MAX_NUMBER_OF_RECORDS);
    info.nsamples = zeros(1, MAX_NUMBER_OF_RECORDS);
    
    if version >= 0.2
        info.recNum = zeros(1, MAX_NUMBER_OF_RECORDS);
    end
    
    current_sample = 0;
    
    RECORD_SIZE = 10 + SAMPLES_PER_RECORD*2 + 10; % size of each continuous record in bytes
    if version >= 0.2
        RECORD_SIZE = RECORD_SIZE + 2; % include recNum
    end
    
    while ftell(fid) + RECORD_SIZE <= filesize % at least one record remains
        
        go_back_to_start_of_loop = 0;
        
        index = index + 1;
        
        if (version >= 0.1)
            timestamp = fread(fid, 1, 'int64', 0, 'l');
            nsamples = fread(fid, 1, 'uint16',0,'l');
            
            
            if version >= 0.2
                recNum = fread(fid, 1, 'uint16');
            end
            
        else
            timestamp = fread(fid, 1, 'uint64', 0, 'l');
            nsamples = fread(fid, 1, 'int16',0,'l');
        end
        
        
        if nsamples ~= SAMPLES_PER_RECORD && version >= 0.1
            
            disp('  Found corrupted record...searching for record marker.');
            
            % switch to searching for record markers
            
            last_ten_bytes = zeros(size(RECORD_MARKER));
            
            for bytenum = 1:RECORD_SIZE*5
                
                byte = fread(fid, 1, 'uint8');
                
                last_ten_bytes = circshift(last_ten_bytes,-1);
                
                last_ten_bytes(10) = double(byte);
                
                if last_ten_bytes(10) == RECORD_MARKER(end)
                    
                    sq_err = sum((last_ten_bytes - RECORD_MARKER).^2);
                    
                    if (sq_err == 0)
                        disp(['   Found a record marker after ' int2str(bytenum) ' bytes!']);
                        go_back_to_start_of_loop = 1;
                        break; % from 'for' loop
                    end
                end
            end
            
            % if we made it through the approximate length of 5 records without
            % finding a marker, abandon ship.
            if bytenum == RECORD_SIZE*5
                
                disp(['Loading failed at block number ' int2str(index) '. Found ' ...
                    int2str(nsamples) ' samples.'])
                
                break; % from 'while' loop
                
            end
            
            
        end
        
        if ~go_back_to_start_of_loop
            
            block = fread(fid, nsamples, 'int16', 0, 'b'); % read in data
            
            fread(fid, 10, 'char*1'); % read in record marker and discard
            
            data(current_sample+1:current_sample+nsamples) = block;
            
            current_sample = current_sample + nsamples;
            
            info.ts(index) = timestamp;
            info.nsamples(index) = nsamples;
            
            if version >= 0.2
                info.recNum(index) = recNum;
            end
            
        end
        
    end
    
    % crop data to the correct size
    data = data(1:current_sample);
    info.ts = info.ts(1:index);
    info.nsamples = info.nsamples(1:index);
    
    if version >= 0.2
        info.recNum = info.recNum(1:index);
    end
    
    % convert to microvolts
    data = data.*info.header.bitVolts;
    
    timestamps = nan(size(data));
    
    current_sample = 0;
    
    if version >= 0.1
        
        for record = 1:length(info.ts)
            
            ts_interp = info.ts(record):info.ts(record)+info.nsamples(record);
            
            timestamps(current_sample+1:current_sample+info.nsamples(record)) = ts_interp(1:end-1);
            
            current_sample = current_sample + info.nsamples(record);
        end
    else % v0.0; NOTE: the timestamps for the last record will not be interpolated
        
        for record = 1:length(info.ts)-1
            
            ts_interp = linspace(info.ts(record), info.ts(record+1), info.nsamples(record)+1);
            
            timestamps(current_sample+1:current_sample+info.nsamples(record)) = ts_interp(1:end-1);
            
            current_sample = current_sample + info.nsamples(record);
        end
        
    end
    
    
    %-----------------------------------------------------------------------
    %--------------------------- SPIKE DATA --------------------------------
    %-----------------------------------------------------------------------
    
elseif strcmp(filetype, 'spikes')
    
    disp('Loading spikes file...');
    
    index = 0;
    
    hdr = fread(fid, NUM_HEADER_BYTES, 'char*1');
    eval(char(hdr'));
    info.header = header;
    
    if (isfield(info.header, 'version'))
        version = info.header.version;
    else
        version = 0.0;
    end
    
    num_channels = info.header.num_channels;
    num_samples = 40; % **NOT CURRENTLY WRITTEN TO HEADER**
    
    % pre-allocate space for spike data
    data = zeros(MAX_NUMBER_OF_SPIKES, num_samples, num_channels);
    timestamps = zeros(MAX_NUMBER_OF_SPIKES, 1);
    info.source = zeros(MAX_NUMBER_OF_SPIKES, 1);
    info.gain = zeros(MAX_NUMBER_OF_SPIKES, num_channels);
    info.thresh = zeros(MAX_NUMBER_OF_SPIKES, num_channels);
    
    if (version >= 0.4)
        info.sortedId = zeros(MAX_NUMBER_OF_SPIKES, num_channels);
    end
    
    if (version >= 0.2)
        info.recNum = zeros(MAX_NUMBER_OF_SPIKES, 1);
    end
    
    
    current_spike = 0;
    last_percent=0;
    
    while ftell(fid) + 512 < filesize % at least one record remains
        
        current_spike = current_spike + 1;
        
        current_percent= round(100* ((ftell(fid) + 512) / filesize));
        if current_percent >= last_percent+10
            last_percent=current_percent;
            fprintf(' %d%%',current_percent);
        end
        
        idx = 0;
        
        % read in event type (1 byte)
        event_type = fread(fid, 1, 'uint8'); % always equal to 4; ignore
        
        idx = idx + 1;
        
        if (version == 0.3)
            event_size = fread(fid, 1, 'uint32', 0, 'l');
            idx = idx + 4;
            ts = fread(fid, 1, 'int64', 0, 'l');
            idx = idx + 8;
        elseif (version >= 0.4)
            timestamps(current_spike) = fread(fid, 1, 'int64', 0, 'l');
            idx = idx + 8;
            ts_software = fread(fid, 1, 'int64', 0, 'l');
            idx = idx + 8;
        end
        
        if (version < 0.4)
            if (version >= 0.1)
                timestamps(current_spike) = fread(fid, 1, 'int64', 0, 'l');
            else
                timestamps(current_spike) = fread(fid, 1, 'uint64', 0, 'l');
            end
            
            idx = idx + 8;
        end
        
        info.source(current_spike) = fread(fid, 1, 'uint16', 0, 'l');
        
        idx = idx + 2;
        
        num_channels = fread(fid, 1, 'uint16', 0, 'l');
        num_samples = fread(fid, 1, 'uint16', 0, 'l');
        
        idx = idx + 4;
        
        if num_samples < 1 || num_samples > 10000
            disp(['Loading failed at block number ' int2str(current_spike) '. Found ' ...
                int2str(num_samples) ' samples.'])
            break;
        end
        
        if (version >= 0.4)
            info.sortedId(current_spike) = fread(fid, 1, 'uint16', 0, 'l');
            electrodeId = fread(fid, 1, 'uint16', 0, 'l');
            channel = fread(fid, 1, 'uint16', 0, 'l');
            color = fread(fid, 3, 'uint8', 0, 'l');
            pcProj = fread(fid, 2, 'single');
            sampleFreq = fread(fid, 1, 'uint16', 0, 'l');
            idx = idx + 19;
        end
        
        waveforms = fread(fid, num_channels*num_samples, 'uint16', 0, 'l');
        
        idx = idx + num_channels*num_samples*2;
        
        wv = reshape(waveforms, num_samples, num_channels);
        
        if (version < 0.4)
            channel_gains = fread(fid, num_channels, 'uint16', 0, 'l');
            idx = idx + num_channels * 2;
        else
            channel_gains = fread(fid, num_channels, 'single');
            idx = idx + num_channels * 4;
        end
        
        info.gain(current_spike,:) = channel_gains;
        
        channel_thresholds = fread(fid, num_channels, 'uint16', 0, 'l');
        idx = idx + num_channels * 2;
        
        info.thresh(current_spike,:) = channel_thresholds;
        
        if version >= 0.2
            info.recNum(current_spike) = fread(fid, 1, 'uint16', 0, 'l');
            idx = idx + 2;
        end
        
        data(current_spike, :, :) = wv;
        
    end
    fprintf('\n')
    for ch = 1:num_channels % scale the waveforms
        data(:, :, ch) = double(data(:, :, ch)-32768)./(channel_gains(ch)/1000);
    end
    
    data = data(1:current_spike,:,:);
    timestamps = timestamps(1:current_spike);
    info.source = info.source(1:current_spike);
    info.gain = info.gain(1:current_spike);
    info.thresh = info.thresh(1:current_spike);
    
    if version >= 0.2
        info.recNum = info.recNum(1:current_spike);
    end
    
    if version >= 0.4
        info.sortedId = info.sortedId(1:current_spike);
    end
    
else
    
    error('File extension not recognized. Please use a ''.continuous'', ''.spikes'', or ''.events'' file.');
    
end

fclose(fid); % close the file

if (isfield(info.header,'sampleRate'))
    if ~ischar(info.header.sampleRate)
        timestamps = timestamps./info.header.sampleRate; % convert to seconds
    end
end

end


function filesize = getfilesize(fid)

fseek(fid,0,'eof');
filesize = ftell(fid);
fseek(fid,0,'bof');

end