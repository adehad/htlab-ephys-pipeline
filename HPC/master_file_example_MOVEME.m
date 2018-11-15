% default options are in parenthesis after the comment

addpath(genpath('KiloSort/')) % path to kilosort folder
addpath(genpath('npy-matlab/')) % path to npy-matlab scripts

pathToYourConfigFile = ''; % take from Github folder and put it somewhere else (together with the master_file)
run(fullfile(pathToYourConfigFile, 'StandardConfig_MOVEME.m'))

tic; % start timer
%
if ops.GPU     
    gpuDevice(1); % initialize GPU (will erase any existing GPU arrays)
end

if strcmp(ops.datatype , 'openEphys')
   ops = convertOpenEphysToRawBInary(ops);  % convert data, only for OpenEphys
end
%
[rez, DATA, uproj] = preprocessData(ops); % preprocess data and extract spikes for initialization
rez                = fitTemplates(rez, DATA, uproj);  % fit templates iteratively
rez                = fullMPMU(rez, DATA);% extract final spike times (overlapping extraction)

% AutoMerge. rez2Phy will use for clusters the new 5th column of st3 if you run this)
%     rez = merge_posthoc2(rez);

% save matlab results file
save(fullfile(ops.root,  'rez.mat'), 'rez', '-v7.3');

% save python results file for Phy
mkdir preAutoMerge
rezToPhy(rez, [ops.root,'preAutoMerge/']);
% To prevent conflicts that arise due to the .phy folder that is
% created in the same location as the binary file
if ismac
    % Code to run on Mac platform
    fprintf(['Mac currently not supported. To use Phy: \n', ...
    'Please move the binary file to the folder \n', ...
    'that you want to visualise \n'])
elseif isunix
    % Code to run on Linux platform
    command = ['ln',' ',ops.fbinary,' ','preAutoMerge/',ops.fbinary];
    [status,cmdout] = system(command)
elseif ispc
    % Code to run on Windows platform
    command = ['mklink /H',' ','preAutoMerge\',ops.fbinary,' ',ops.fbinary];
    [status,cmdout] = system(command)
else
    disp('Platform not supported')
end


%%
% AUTO MERGES 
% after spending quite some time with Phy checking on the results and understanding the merge and split functions, 
% come back here and run Kilosort's automated merging strategy. This block
% will overwrite the previous results and python files. Load the results in
% Phy again: there should be no merges left to do (with the default simulation), but perhaps a few splits
% / cleanup. 
% Kilosort's AUTO merges should not be confused with the "best" merges 

rez = merge_posthoc2(rez);

% save python results file for Phy
mkdir postAutoMerge
rezToPhy(rez, [ops.root,'postAutoMerge/']);
% To prevent conflicts that arise due to the .phy folder that is
% created in the same location as the binary file
if ismac
    % Code to run on Mac platform
    fprintf(['Mac currently not supported. To use Phy: \n', ...
    'Please move the binary file to the folder \n', ...
    'that you want to visualise \n'])
elseif isunix
    % Code to run on Linux platform
    command = ['ln',' ',ops.fbinary,' ','postAutoMerge/',ops.fbinary];
    [status,cmdout] = system(command)
elseif ispc
    % Code to run on Windows platform
    command = ['mklink /H',' ','postAutoMerge\',ops.fbinary,' ',ops.fbinary];
    [status,cmdout] = system(command)
else
    disp('Platform not supported')
end





% remove temporary file
delete(ops.fproc);

display('Done')
