%% Wrapper File for Extracting Information from Binary files & Sorting Program Outputs
% Use Ctrl+Enter (Windows) to run section by section
%% Add Functions to Path
addpath(genpath('C:\Users\Daniel\Documents\GitHub\htlab-ephys-pipeline\MATLAB-Analysis\spikeSorting')) % Path to spikeSorting folder of git
%% SECTION 1: Extracts spike times to _sortedmat
sortedType = 'kilosort';            
sortOutputFolder = 'C:\PATH\TO\THE\SORTED\FOLDER\181108\preAutoMerge';
file = '181108';
startTrial = 1;
mergedInfoCSV = [];
clusterType = "unsorted"; % good, unsorted, noise or [] (empty) (when there has been no manual sorting, there's only unsorted). 
                          % Arrays accepted e.g. ["good";"unsorted";"MUA"]
filename.sortOutputAll = extractTrialUnits(sortedType,...           % sorting Program used
                                       sortOutputFolder, ...    % location of sorting output
                                       file, ....               % ['YYMMDD' ; 'XYZ.kwik'] - XYZ.kwik is only if you are using klusta
                                       startTrial, ...          % YYMMDD_X - where is is starting number
                                       mergedInfoCSV, ...       % Path + Name of merge_info.csv
                                       clusterType);            % ['XYZ' ; 'ABC'] - rows containing different cluster types to keep - note: depends on your manual clustering - leave empty if you want to keep all of them
% Select a single sortOutput
filename.sortOutput = [sortOutputFolder, filesep, filename.sortOutputAll(1,:)];   
%% SECTION 2: Establish Metafile struct
filename.folder = 'C:\Users\Daniel\Box Sync\DragonVision_DanielKo\Data\Ephys\181116\'; 
filename.binary=[filename.folder, filesep, '2018-11-16_16-23-50\2018-11-16_16-23-50_padded.bin'];
filename.sortOutput=[filename.folder, filesep, '181116_05_sorted.mat'];
if isempty(filename.binary)
    [fileName, filePath] = uigetfile('*.bin','Select experiment raw binary data .bin file:');
    filename.binary = [filePath filesep fileName];
end
if isempty(filename.sortOutput)
    [fileName,filePath] = uigetfile('*.mat','Select experiment _sorted.mat file:');
    filename.sortOutput = [filePath filesep fileName];
end

%%% IF YOU HAVE A .meta FILE (telemetry)
    %[m, fpath, mfile] = readMetafile2('150526__MovingObjects_1.meta','C:\PATH\TO\THE\METAFILE\150526\Tetrode test data\');
    %[m, fpath, mfile] = readMetafile(); % GUI file picker
    %m.metafile = mfile;
    %m.metapath = fpath;
%%%

%%% IF YOU DO NOT HAVE A .meta FILE
     %m.nChans was here
     m.sRateHz = 30e3;   % sampling frequency
%%%

m.dbytes    = 2; % byte size of data - i.e. int16 is 2 bytes
m.msec      = m.sRateHz/1000; % conversion factor from ms time to sample number

%% SECTION 3: Extract Waveforms
m.nChans    = 4;            % number of channels
m.ech       = 1:m.nChans-1; % ephys channel(s) is everything except the last

[m,s] = extractTrialUnitWaves(filename.binary, ... % Binary File
                      filename.sortOutput, ...  % _sorted.mat file
                      m, ...                    % metafile struct, m
                      0, ...                    % 1: if you want to do secondary template matching
                      []);              % filename to store output, leave as [] if you don't want to save
                                                           
%% SECTION 4: Saves Output for a given filename.sortOutput
% renamed YYMMDD_XX_sorted.mat -> YYMMDD_XX_extracted.mat
save([filename.sortOutput(1:end-11),'_extracted.mat'], 'm', 's', 'filename')
