%% Wrapper File for Extracting Information from Binary files & Sorting Program Outputs
%% Add Functions to Path
addpath(genpath('C:\Users\Daniel\Documents\GitHub\htlab-ephys-pipeline\MATLAB-Analysis\visualAnalysis'))
%% SECTION 1: Extract PD data
m.nChans    = 1;        % number of channels
m.pdch      = 1;        %assume pd is last ch
m.fps       = 180;      % (projector frame rate)*3  (*3 if B&W)
% m.pdthr = 3e3; % Can set now, or comment out to set graphically in function (telemetry)
m.pdthr     = 3;        % OpenEphys projection PD 

m = extractTrialADC_PD(filename.rawADC, ... % Binary File
                        m, ....     % metafile struct, m
                        [] ); % filename to store output, leave as [] if you don't want to save

%% SECTION 2: Add stimulus and projection data
nLoops = 5;
pdDiffThreshold = 1.5e3;

stim = addProj2M(m, nLoops, pdDiffThreshold);

%% SECTION 2.5: Save .mat
save('181116_05_visual.mat', 'm', 's', 'stim', 'filename')

%% SECTION 2.9999: Options
opt.drawingMode
opt.units
opt.
%% SECTION 3: Raster plot