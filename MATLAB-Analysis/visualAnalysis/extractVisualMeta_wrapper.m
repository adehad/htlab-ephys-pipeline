%% Wrapper File for Extracting Information from Binary files & Extracting Projection Settings
%% Add Functions to Path
addpath(genpath('E:\GitHub\htlab-ephys-pipeline\MATLAB-Analysis\visualAnalysis'))
%addpath(genpath('C:\Users\Daniel\Documents\GitHub\htlab-ephys-pipeline\MATLAB-Analysis\visualAnalysis'))

%% SECTION 0: Load spike data .mat - only if you haven't loaded it yet
uiopen('*_extracted.mat');

%% SECTION 1: Extract PD data
% path to photodiode data _ADC.bin file
filename.rawADC=[];%'K:\Ephys\181207\1\2018-12-07_19-26-39\2018-11-16_16-23-50_ADC.bin';
if isempty(filename.rawADC)
    [fileName, filePath] = uigetfile('*.bin','Select experiment _ADC.bin file:');
    filename.rawADC = [filePath fileName];
end

m.nChans    = 1;        % number of channel in adc bin
m.pdch      = 1;        % the index of photodiode channel
m.fps       = 180;      % projector subframerate

% set these pd event thresholds now or set graphically in extractTrialADC_PD
% m.pdthr = 3e3; % telemetry
m.pdthr     = 3;        % OpenEphys projection PD 

m = extractTrialADC_PD(filename.rawADC, ... % Binary File
                        m, ....     % metafile struct, m
                        []); % filename to store output, leave as [] if you don't want to save

%% SECTION 2: Add stimulus and projection data
% plot interval between pd events so you can choose pdDIffThreshold
plot(diff(m.pd)); ylim([0 2e5]);

%% SECTION 2: Add stimulus and projection data
clear stim
%%%%%%%%%%%%%%%%%%%%%%% IMPORTANT
stim.StimGL_nloops = 16;     % number of loops in stimulus
pdDiffThreshold = 1500;    % threshold of interval between loops
%%%%%%%%%%%%%%%%%%%%%%%

stim.C=0;
stim.D=13.5;                % y Distance of DF head from centre of screen in mm
stim.E=22.5;                % z Distance of the DF head out of the screen in mm
stim.W=67;                  % Width of the projection in mm
stim.H=91;                  % Height of the projection in mm

stim.xPix = 480;            % xaxis pixel count
stim.yPix = 640;            % yaxis pixel count

stim = addProj2M(m, stim, pdDiffThreshold);

%% SECTION 3: Save .mat
save('181017_06_visual_pruned.mat', 'm', 's', 'stim', 'filename')