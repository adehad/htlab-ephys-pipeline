%%%% Split_trials
%%%% A. Haddad, H. Lin
%%%% INPUT: spike_time.npy, spike_clusters.npy
%%%% This script takes the npy files from running KiloSort on merged experiment data
%%%% and split data back to individual trials with different clusters
%%%% 
%%% Required External Functions:
%%% readOurCSV  readNPY readNPYheader

%% Check Required Functions
addpath(genpath('requiredFunctions')) % path to folder with functions
reqFuncStr = '';
if ~exist('readNPY.m','file');   reqFuncStr=strcat(reqFuncStr,'readNPY.m (npy-matlab github), ');end4if ~exist('readNPY.m','file');   reqFuncStr=strcat(reqFuncStr,'readNPY.m (npy-matlab github), ');end
if ~exist('readNPYheader.m','file');   reqFuncStr=strcat(reqFuncStr,'readNPYheader.m (npy-matlab github), ');end
if ~exist('readOurCSV.m','file');    reqFuncStr=strcat(reqFuncStr,'readOurCSV.m, '); end
% if ~exist('getWaveForms.m','file');        reqFuncStr=strcat(reqFuncStr,'getWaveForms.m (spikes github), ');      end
% if ~exist('rasterplot.m','file');        reqFuncStr=strcat(reqFuncStr,'rasterplot.m, ');      end

if ~strcmp(reqFuncStr,'')
error(['The following functions could not be found in the search path:', newline, reqFuncStr])
end
%% trial info
% file_dir = 
files = ['131210_01'; '131210_02'];
startTrial = 4; % e.g. if merged trials 4 to 6, startTrial = 4

%% parse clusters
% This part of the script must be run insde the kilosort output directory
% need to set it to the correct directory first
% addpath(genpath('C:\Users\Adehad\Desktop\dragonFlyUROP\Code\KiloSort\npy-matlab-master')) % path to npy-matlab scripts

spike_times = readNPY('spike_times.npy');
spike_clusters = readNPY('spike_clusters.npy');
clusters=sort(unique(spike_clusters));

% for cc=1:length(clusters) % note that cluster 0 is always leftovers
%      eval(['unit_',num2str(cc-1),'= spike_times(find(spike_clusters == clusters(cc)));']) 
% 
% end

%% split back into trials
% load('YYMMDD_merge_info.txt') % still need to import it properly
fileDetails = readOurCSV(['131210_merge46_info','.csv']);
fileDetails.filename % a list of all the file names
fileDetails.index % a list of where in the merged dataset each file *ends*
fileDetails.samples % a list where each element is the length in samples of each file

split_point = [0; fileDetails.index]; 


n_trial=length(split_point);
%%
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
        eval(['s.unit_',num2str(clusters(cc)),'= trial_spikes(trial_clusters == clusters(cc));'])
    end
    
    if tt>9    % Create 2 digit IDss
        trialID = num2str(startTrial+tt-1);
    else
        trialID = ['0',num2str(startTrial+tt-1)];
    end    
    save([files(1,1:6),'_',trialID,'_sorted.mat'],'s','trial_spikes','trial_clusters') % save structure s, trial_cluster, and trial_spikes to each individiaul trial as yymmdd_nn_sorted.mat
    disp(['saved data for trial: ',trialID])
end

