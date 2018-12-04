%% Wrapper File for Extracting Information from Binary files & Sorting Program Outputs
%% Add Functions to Path
addpath(genpath('C:\Users\Daniel\Documents\GitHub\htlab-ephys-pipeline\MATLAB-Analysis\spikeSorting'))
%% SECTION 1: Extracts spike times to _sortedmat
sortedType = 'kilosort';
sortOutputFolder = 'C:\PATH\TO\THE\SORTED\FOLDER\181108\preAutoMerge';
file = '181108';
startTrial = 1;
csvName = [];
clusterType = 'unsorted'; %good, unsorted, noise or empty (when there has been no manual sorting, theres only unsorted). Arrays accepted
filename.sortOutput= extractTrialUnits(sortedType,...           % sorting Program used
                                       sortOutputFolder, ...    % location of sorting output
                                       file, ....               % ['YYMMDD' ; 'XYZ.kwik'] - XYZ.kwik is only if you are using klusta
                                       startTrial, ...          % YYMMDD_X - where is is starting number
                                       csvName, ...             % Name of merge_info csv
                                       clusterType);            % ['XYZ' ; 'ABC'] - rows containing different cluster types to keep - note: depends on your manual clustering - leave empty if you want to keep all of them

filename.sortOutput = [sortOutputFolder, filename.sortOutput];   
%% SECTION 2: Establish Metafile struct
filename.raw=['C:\Users\Daniel\Box Sync\DragonVision_DanielKo\Data\Ephys\181116\2018-11-16_16-23-50\2018-11-16_16-23-50_padded.bin'];
filename.rawADC=['C:\Users\Daniel\Box Sync\DragonVision_DanielKo\Data\Ephys\181116\2018-11-16_16-23-50\2018-11-16_16-23-50_ADC.bin'];
filename.sortOutput=['C:\Users\Daniel\Box Sync\DragonVision_DanielKo\Data\Ephys\181116\181116_05_sorted.mat'];

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

[m,s] = extractTrialUnitWaves(filename.raw, ... % Binary File
                      filename.sortOutput, ...  % _sorted.mat file
                      m, ...                    % metafile struct, m
                      0, ...                    % 1: if you want to do secondary template matching
                      []);              % filename to store output, leave as [] if you don't want to save
                                                           
%% SECTION 4: Extract PD data
m.nChans    = 1;        % number of channels
m.pdch      = 1;        %assume pd is last ch
m.fps       = 180;      % (projector frame rate)*3  (*3 if B&W)
% m.pdthr = 3e3; % Can set now, or comment out to set graphically in function (telemetry)
m.pdthr     = 3;        % OpenEphys projection PD 

m = extractTrialADC_PD(filename.rawADC, ... % Binary File
                        m, ....     % metafile struct, m
                        [] ); % filename to store output, leave as [] if you don't want to save

%%
save('181116_05.mat', 'm', 's', 'filename')
