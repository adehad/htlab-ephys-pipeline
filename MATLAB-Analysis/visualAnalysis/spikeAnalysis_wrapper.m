%% Add functions to path
addpath(genpath('E:\GitHub\htlab-ephys-pipeline\MATLAB-Analysis\visualAnalysis'))
%addpath(genpath('C:\Users\Daniel\Documents\GitHub\htlab-ephys-pipeline\MATLAB-Analysis\visualAnalysis'))

% Better figure exporting
addpath(genpath('K:\MATLAB\altmany-export_fig-26eb699'))
%% SECTION 0: Load spike and pd data
uiopen('*_visual.mat');

%% SECTION 1: Options
opt.preName = 'BRININDAHEAT'; % first part of name of saved figures
opt.tsdnLatency = 0; % latency of tsdns in ms
selectUnits = 4; % array of units to plot {'all [X,Y,Z]}

opt.rasterBinSize = 100; % for raster: bin size of spike rate histogram in ms

opt.drawingMode = 'none'; % for heatmap: how to draw orientation distribution within bin {'none, 'arrows', 'blob'}
opt.units = 'angle'; % for heatmap and topography: units on sides of heatmap {'pixels' 'angle'}
opt.sqSize = 10; % for heatmap and topography: sizes of square bins
opt.halfOffset = 1;
opt.discardCorner = 1; % for heatmap, polar and topography: whether to ignore out of bounds data {1, 0}

%% SECTION 2: Raster plot
plotRaster(m,s,stim,selectUnits,opt,1); % last argument is save? {0 1 2 = export_fig}

%% SECTION 3: Heatmap
plotHeatmap(m,s,stim,selectUnits,opt,0); % last argument is save? {0 1 2 = export_fig}

%% SECTION 4: Polar plot
plotPolar(m,s,stim,selectUnits,opt,1); % last argument is save? {0 1 2 = export_fig}

%% SECTION 5: Orientation topography
plotOrientationTopography(m,s,stim,selectUnits,opt,1); % last argument is save? {0 1 2 = export_fig}

%% SECTION 6: Trajectory Histogram
[n1, q1, b1] = plotTrajHist(m,s,stim,selectUnits,opt,1); % last argument is save? {0 1 2 = export_fig}