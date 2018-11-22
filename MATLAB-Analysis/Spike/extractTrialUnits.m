function extractUnitsFileName = extractTrialUnits(sortedType,sortOutputFolder, file, startTrial, csvName, clusterType)
%%% COPY PASTE JOB FOR NOW
% sortedType = klusta or kilosort
% file = [ YYMMDD ; KWIK FILE]      % KWIK file only needed if 'klusta'
% csvName = merge_info.csv % can also contain the path e.g. C:\folder\merge_info.csv
% clusterType = array of ones to extract e.g. good, unsorted
% Example Usage:


%% Check Required Functions
% Appended to the end of this fuction
% readNPY   readNPYheader   readOurCSV  getKwikData      

%% read the merge_info csv file (csvName can contain
if ~isempty(csvName)        % if name is not empty
    fileDetails = readOurCSV(csvName);
    fileDetails.filename % a list of all the file names
    fileDetails.index % a list of where in the merged dataset each file *ends*
    fileDetails.samples % a list where each element is the length in samples of each file

    split_point = [0; fileDetails.index]; 
else % Assume file is not merged - just a single experiment
    split_point = [0; Inf]; % So we extract all spike_times for when we use non-merged data ?
end

%% cd to the folder with the sorted data
    pwdStore = pwd; % first store where we currently are 
    cd(sortOutputFolder)
    
    clear s % reset the structure
    
%% Cluster extraction based on sorting program
if strcmpi(sortedType,'kilosort') % case insensitive srtcmp
    %% Read KiloSort Files
    % This part of the script must be run insde the kilosort output directory
    spike_times = readNPY('spike_times.npy');
    spike_clusters = readNPY('spike_clusters.npy');
    
    % To overcome KiloSort issue of duplication  
    duplAmount = length(spike_times)/length(unique(spike_times)); %  number of times spike time is duplicated
    spike_times = spike_times(1:duplAmount:end); 
    spike_clusters = spike_clusters(1:duplAmount:end);
    %% Find Cluster groupings
    if isfile('cluster_groups.csv') % i.e. file exists
        cluster_groups = tdfread('cluster_groups.csv');
    else
        error(['cluster_groups.csv does not exist in the provided folder', newline, ...
               sortOutputFolder,newline,'Make sure you Ctrl+S after you cluster in Phy']);
    end
    %% Store info for desiredclusters
    clusters = [];
    if ~isempty(clusterType) % cluster types specified by input arguments - e.g. ['good';'unsorted']
        for cc=1:size(clusterType,1)
            for ii=1:size(cluster_groups.cluster_id,1)
                if contains( cluster_groups.group(ii,:), clusterType(cc,:),'IgnoreCase',true )
                    clusters = [clusters cluster_groups.cluster_id(ii)]; %accumulate only the good clusters
                end
            end
        end
        s.cluster_groups = cluster_groups; % preserving the original allocations by storing this struct
        s.clusters = clusters;      % storage of what clusters were kepts
    else % if no cluster types specified, just take them all
        clusters=sort(unique(spike_clusters));
        if ~isempty(cluster_groups)
            s.cluster_groups = cluster_groups;
        end
    end


elseif strcmpi(sortedType,'klusta') % case insensitive srtcmp
    %% Read Kwik files to obtain cluster info
    kwikFileName = file(2,:); % second row of file must be the kwik file name
    kwik_d=getKwikData(kwikFileName);
    spike_times = hdf5read(kwikFileName, '/channel_groups/0/spikes/time_samples');
    spike_clusters = hdf5read(kwikFileName, '/channel_groups/0/spikes/clusters/main');
    clusters = [];
    for cc=1:size(clusterType,1)
        for ii=1:length(kwik_d)
            clu_type=cell2mat(kwik_d(ii).id(1,1));
            if strcmpi(clu_type, clusterType(cc,:)) % case insensitive compare
                clusters = [clusters kwik_d(ii).icell]; %accumulate only the good clusters
            end
        end
    end
    
    s.clusters = clusters;      % storage of what clusters were kepts
else
    error(['sortedType provided does not match the accepted :', newline, ...
              ' ''kilosort'' ', newline,' ''klusta'' ' ]);
end

%% Using trial information, Associate each spike time with each trial
% trials refer to the experiment binary files that were merged together

n_trial=length(split_point);

% split back into trials
for tt=1:n_trial-1 % go through each trial
    trial_spikes = spike_times(spike_times>split_point(tt) & spike_times <= split_point(tt+1));
    trial_spikes = trial_spikes - split_point(tt);
            % spike_times within that trial (i.e. between split times)
            % We need to index them as if they were individual files, hence
            % the subtraction of the other index
    trial_clusters = spike_clusters(spike_times>split_point(tt) & spike_times <= split_point(tt+1));
        % corresponding cluster assignments within that trial      
    for cc=1:length(clusters) % go through each cluster
        s.(['unit_',num2str(clusters(cc),'%02i')]) =  trial_spikes(trial_clusters == clusters(cc));
        s.units{clusters(cc)} = trial_spikes(trial_clusters == clusters(cc));
    end

    % Create 2 digit IDss
        trialID = num2str(startTrial+tt-1,'%02i');
   
    extractUnitsFileName = [file(1,1:6),'_',trialID,'_sorted.mat'];
    save(extractUnitsFileName,'s','trial_spikes','trial_clusters') % save structure s, trial_cluster, and trial_spikes to each individiaul trial as yymmdd_nn_sorted.mat
    disp(['saved data for trial: ',trialID]);
    
    extractUnitsFileName = [file(1,1:6),'_',trialID,'_sorted.mat'];
    
end
    

cd(pwdStore) % return to original folder
end


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


function [arrayShape, dataType, fortranOrder, littleEndian, totalHeaderLength, npyVersion] = readNPYheader(filename)
% function [arrayShape, dataType, fortranOrder, littleEndian, ...
%       totalHeaderLength, npyVersion] = readNPYheader(filename)
%
% parse the header of a .npy file and return all the info contained
% therein.
%
% Based on spec at http://docs.scipy.org/doc/numpy-dev/neps/npy-format.html

fid = fopen(filename);

% verify that the file exists
if (fid == -1)
    if ~isempty(dir(filename))
        error('Permission denied: %s', filename);
    else
        error('File not found: %s', filename);
    end
end

try
    
    dtypesMatlab = {'uint8','uint16','uint32','uint64','int8','int16','int32','int64','single','double', 'logical'};
    dtypesNPY = {'u1', 'u2', 'u4', 'u8', 'i1', 'i2', 'i4', 'i8', 'f4', 'f8', 'b1'};
    
    
    magicString = fread(fid, [1 6], 'uint8=>uint8');
    
    if ~all(magicString == [147,78,85,77,80,89])
        error('readNPY:NotNUMPYFile', 'Error: This file does not appear to be NUMPY format based on the header.');
    end
    
    majorVersion = fread(fid, [1 1], 'uint8=>uint8');
    minorVersion = fread(fid, [1 1], 'uint8=>uint8');
    
    npyVersion = [majorVersion minorVersion];
    
    headerLength = fread(fid, [1 1], 'uint16=>uint16');
    
    totalHeaderLength = 10+headerLength;
    
    arrayFormat = fread(fid, [1 headerLength], 'char=>char');
    
    % to interpret the array format info, we make some fairly strict
    % assumptions about its format...
    
    r = regexp(arrayFormat, '''descr''\s*:\s*''(.*?)''', 'tokens');
    dtNPY = r{1}{1};    
    
    littleEndian = ~strcmp(dtNPY(1), '>');
    
    dataType = dtypesMatlab{strcmp(dtNPY(2:3), dtypesNPY)};
        
    r = regexp(arrayFormat, '''fortran_order''\s*:\s*(\w+)', 'tokens');
    fortranOrder = strcmp(r{1}{1}, 'True');
    
    r = regexp(arrayFormat, '''shape''\s*:\s*\((.*?)\)', 'tokens');
    shapeStr = r{1}{1}; 
    arrayShape = str2num(shapeStr(shapeStr~='L'));

    
    fclose(fid);
    
catch me
    fclose(fid);
    rethrow(me);
end


end


function data = readNPY(filename)
% Function to read NPY files into matlab. 
% *** Only reads a subset of all possible NPY files, specifically N-D arrays of certain data types. 
% See https://github.com/kwikteam/npy-matlab/blob/master/npy.ipynb for
% more. 
%

[shape, dataType, fortranOrder, littleEndian, totalHeaderLength, ~] = readNPYheader(filename);

if littleEndian
    fid = fopen(filename, 'r', 'l');
else
    fid = fopen(filename, 'r', 'b');
end

try
    
    [~] = fread(fid, totalHeaderLength, 'uint8');

    % read the data
    data = fread(fid, prod(shape), [dataType '=>' dataType]);
    
    if length(shape)>1 && ~fortranOrder
        data = reshape(data, shape(end:-1:1));
        data = permute(data, [length(shape):-1:1]);
    elseif length(shape)>1
        data = reshape(data, shape);
    end
    
    fclose(fid);
    
catch me
    fclose(fid);
    rethrow(me);
end


end


function chans = getKwikData(kwikFile)
igroup=0; basename = 'XXX';
spkTimes = hdf5read(kwikFile, ['/channel_groups/' num2str(igroup) '/spikes/time_samples']);
spkClus = hdf5read(kwikFile, ['/channel_groups/' num2str(igroup) '/spikes/clusters/main']);
cellIDs = unique(spkClus);
ncells  = length(cellIDs);

% spkClus(spkTimes>endTime) = [];
% spkTimes(spkTimes>endTime) = [];
% spkTimes = double(spkTimes);
% spkTimes = spkTimes - startTime;
% spkClus(spkTimes<0) = [];
% spkTimes(spkTimes<0) = [];

temp =  h5info(kwikFile,  '/recordings/0/');
sampleRate = double((temp.Attributes(3).Value));
spkTimes = double(spkTimes)./sampleRate;
noise_list = [];
for icell = 1:ncells
    chans(icell).spiketimes = spkTimes(spkClus==cellIDs(icell));
%     chans(icell).ichan = igroup;
%     chans(icell).iexp = iexp;
    chans(icell).icell = cellIDs(icell);
    chans(icell).sampleRate = sampleRate;
    
    temp = h5info(kwikFile, ['/channel_groups/' num2str(igroup) '/clusters/main/' num2str(cellIDs(icell)) '/']);
    
    for idx = 1:length(temp.Attributes)
        if strcmp(temp.Attributes(idx).Name, 'cluster_group')
            ilabel = temp.Attributes(idx).Value;
            break
        end
        idx = idx + 1;
    end
    if idx <= length(temp.Attributes)
        ilabel = temp.Attributes(idx).Value;
        temp = h5info(kwikFile, ['/channel_groups/' num2str(igroup) '/cluster_groups/main']);
        labelType = eval(['temp.Groups(' num2str(ilabel+1) ').Attributes(end).Value']);
        if strcmp(labelType,'Noise')
            noise_list = [noise_list icell];
        end
        chans(icell).id = [labelType '_cluster' num2str(cellIDs(icell))];
    else
        chans(icell).id = ['none_cluster' num2str(cellIDs(icell))];
    end
end
chans(noise_list) = [];


end