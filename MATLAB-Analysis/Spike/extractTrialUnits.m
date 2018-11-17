function extractUnitsFileName = extractTrialUnits(sortedType,outputFolder, file, startTrial, csvName, clusterType)
%%% COPY PASTE JOB FOR NOW
% sortedType = klusta or kilosort
% file = list of files ? maybe not needed
%       % MUST BE YYMMDD_XX format
%       % USE KWIK FILE NAME IF 'klusta'
% csvName = merge_info.csv 
% clusterType = array of ones to extract e.g. good, unsorted
% Example Usage:


%%% Required External Functions:
%%% readOurCSV  readNPY readNPYheader
extractUnitsFileName = '';
startTrial = 1;

%% Check Required Functions - Append all to bottom of the file instead ?
addpath(genpath('requiredFunctions')) % path to folder with functions
reqFuncStr = '';
if ~exist('readNPY.m','file');          reqFuncStr=strcat(reqFuncStr,'readNPY.m (npy-matlab github), ');end
if ~exist('readNPYheader.m','file');    reqFuncStr=strcat(reqFuncStr,'readNPYheader.m (npy-matlab github), ');end
if ~exist('readOurCSV.m','file');       reqFuncStr=strcat(reqFuncStr,'readOurCSV.m, '); end
if ~exist('getKwikData.m','file');      reqFuncStr=strcat(reqFuncStr,'getKwikData.m, '); end

if ~strcmp(reqFuncStr,'')
error(['The following functions could not be found in the search path:', newline, reqFuncStr])
end

if ~isempty(csvName)        % if name is not empty
    fileDetails = readOurCSV([csvName,'.csv']);
    fileDetails.filename % a list of all the file names
    fileDetails.index % a list of where in the merged dataset each file *ends*
    fileDetails.samples % a list where each element is the length in samples of each file

    split_point = [0; fileDetails.index]; 
    
elseif ~isempty(file)
    
else % Assume file is not merged - just a single experiment
    split_point = [0; Inf]; % So we extract all spike_times for when we use non-merged data ?
end

if strcmp(sortedType,'kilosort')
    %% parse clusters
    % This part of the script must be run insde the kilosort output directory
    % need to set it to the correct directory first

    spike_times = readNPY('spike_times.npy');
    spike_clusters = readNPY('spike_clusters.npy');
    spike_times = unique(spike_times); % To overcome KiloSort issue of duplication
    
    % Find Cluster groupings
    if isfile('cluster_groups.csv') % i.e. file exists
        cluster_groups = tdfread('cluster_groups.csv');
    end
    clusters = [];
    if ~isempty(clusterType)
        for cc=1:size(clusterType,1)
            for ii=1:size(cluster_groups.cluster_id,1)
                if contains( cluster_groups.group(ii,:), clusterType(cc,:),'IgnoreCase',true )
                    clusters = [clusters cluster_groups.cluster_id(ii)]; %accumulate only the good clusters
                end
            end
        end
        s.cluster_groups = cluster_groups;
    else % if no cluster types specified, just take them all
        clusters=sort(unique(spike_clusters));
        if ~isempty(cluster_groups)
            s.cluster_groups = cluster_groups;
        end
    end


elseif strcmp(sortedType,'klusta')
    %% extract spike-sorting data

    kwikFileName = file;
    kwik_d=getKwikData(kwikFileName);
    spike_times = hdf5read(kwikFileName, '/channel_groups/0/spikes/time_samples');
    spike_clusters = hdf5read(kwikFileName, '/channel_groups/0/spikes/clusters/main');
    clusters = [];
    for ii=1:length(kwik_d)
        clu_type=cell2mat(kwik_d(ii).id(1,1));
        if ~isempty(strfind(clu_type,'Good'))
            clusters = [clusters kwik_d(ii).icell]; %accumulate only the good clusters
        end
    end
else
    error(['sortedType provided does not match the accepted :', newline, ...
              ' ''kilosort'' ', newline,' ''klusta'' ' ]);
end

%% extract trial information

n_trial=length(split_point);

%% split back into trials
for tt=1:n_trial-1 % go through each trial
    trial_spikes = spike_times(spike_times>split_point(tt) & spike_times <= split_point(tt+1));
    trial_spikes = trial_spikes - split_point(tt);
            % spike_times within that trial (i.e. between split times)
            % We need to index them as if they were individual files, hence
            % the subtraction of the other index
    trial_clusters = spike_clusters(spike_times>split_point(tt) & spike_times <= split_point(tt+1));
        % corresponding cluster assignments within that trial     
    clear s % reset the structure 
    for cc=1:length(clusters) % go through each cluster
        s.(['unit_',num2str(cc-1)]) =  trial_spikes(trial_clusters == clusters(cc));
    end

    % Create 2 digit IDss
        trialID = num2str(startTrial+tt-1,'%02i');
   
    save([file(1,1:6),'_',trialID,'_sorted.mat'],'s','trial_spikes','trial_clusters') % save structure s, trial_cluster, and trial_spikes to each individiaul trial as yymmdd_nn_sorted.mat
    disp(['saved data for trial: ',trialID])
end
    


end