%%%% Split_trials2
%%%% A. Haddad, H. Lin
%%%% INPUT: (spike_time.npy, spike_clusters.npy) or (xxxx.kwik) and merge_info.csv
%%%% This script takes the spike-sorted merged experiment data
%%%% and split data back to individual trials with different clusters
%%%% The input can be .npy files from KiloSort or .kwik file from Klusta 
%%% Required External Functions:
%%% readOurCSV readNPY readNPYheader

% function Split_trials2(file_type, files, filename)
% %% Check Required Functions
% addpath(genpath('requiredFunctions')) % path to folder with functions
% reqFuncStr = '';
% if ~exist('readNPY.m','file');   reqFuncStr=strcat(reqFuncStr,'readNPY.m (npy-matlab github), ');end4if ~exist('readNPY.m','file');   reqFuncStr=strcat(reqFuncStr,'readNPY.m (npy-matlab github), ');end
% if ~exist('readNPYheader.m','file');   reqFuncStr=strcat(reqFuncStr,'readNPYheader.m (npy-matlab github), ');end
% if ~exist('readOurCSV.m','file');    reqFuncStr=strcat(reqFuncStr,'readOurCSV.m, '); end
% % if ~exist('getWaveForms.m','file');        reqFuncStr=strcat(reqFuncStr,'getWaveForms.m (spikes github), ');      end
% % if ~exist('rasterplot.m','file');        reqFuncStr=strcat(reqFuncStr,'rasterplot.m, ');      end
% 
% if ~strcmp(reqFuncStr,'')
% error(['The following functions could not be found in the search path:', newline, reqFuncStr])
% end
%% trial info
file_type = 0; 
files = ['131230_01'; '131230_02'; '131230_03'; '131230_04'; '131230_05'];

% %% read data and parse clusters
% if file_type==2
%     spike_times = readNPY('spike_times.npy');
%     spike_clusters = readNPY('spike_clusters.npy');
% elseif file_type==0
%     kwik_d=getKwikData('131230_merged.kwik');
%     spike_times = hdf5read('131230_merged.kwik', '/channel_groups/0/spikes/time_samples');
%     spike_clusters = hdf5read('131230_merged.kwik', '/channel_groups/0/spikes/clusters/main');
% 
% % else
% %     spike_times = hdf5read('131230_merged.kwik', '/channel_groups/0/spikes/time_samples');
% %     spike_clusters = cast(importfile('131230_merged.clu.0'),'int32');
% end

%% extract spike-sorting data
kwikFileName = '131230_merged.kwik';
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
% clusters=sort(unique(spike_clusters));

%% extract trial information
% load('YYMMDD_merge_info.txt') % still need to import it properly
fileDetails = readOurCSV([kwikFileName(1:6),'_merge_info','.csv']);
fileDetails.filename % a list of all the file names
fileDetails.index % a list of where in the merged dataset each file *ends*
fileDetails.samples % a list where each element is the length in samples of each file
split_point = [0; fileDetails.index]; 
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
        eval(['s.unit_',num2str(cc-1),'= trial_spikes(trial_clusters == clusters(cc));'])
    end
    
    if tt>9    % Create 2 digit IDss
        trialID = num2str(tt);
    else
        trialID = ['0',num2str(tt)];
    end    
    save([files(1,1:6),'_',trialID,'_sorted.mat'],'s','trial_spikes','trial_clusters') % save structure s, trial_cluster, and trial_spikes to each individiaul trial as yymmdd_nn_sorted.mat
    disp(['saved data for trial: ',trialID])
end

