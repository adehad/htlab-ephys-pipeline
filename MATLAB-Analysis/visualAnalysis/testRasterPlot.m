%% load s and m
addpath(genpath('C:\Users\Daniel\Documents\GitHub\htlab-ephys-pipeline\MATLAB-Analysis\Spike'))

%% Establish Metafile struct
filename.raw=['C:\Users\Daniel\Box Sync\MEGA_BINARY\181017_06.bin'];
filename.sortOutput=['C:\Users\Daniel\Box Sync\MEGA_BINARY\KiloSort\181017_06_sorted.mat'];

%%% IF YOU HAVE A .meta FILE (telemetry)
    %[m, fpath, mfile] = readMetafile2('150526__MovingObjects_1.meta','C:\PATH\TO\THE\METAFILE\150526\Tetrode test data\');
    %[m, fpath, mfile] = readMetafile(); % GUI file picker
    %m.metafile = mfile;
    %m.metapath = fpath;
%%%

%%% IF YOU DO NOT HAVE A .meta FILE
     m.nChans = 2;       % number of channels
     m.sRateHz = 30e3;   % sampling frequency
%%%
m.pdch      = m.nChans; %assume pd is last ch
m.ech       = 1:m.nChans-1; % ephys channel(s) is everything except the last
m.dbytes    = 2; % byte size of data - i.e. int16 is 2 bytes
m.msec      = m.sRateHz/1000; % conversion factor from ms time to sample number

%% Extract Waveforms
[m,s] = extractTrialUnitWaves(filename.raw, ... % Binary File
                      filename.sortOutput, ...  % _sorted.mat file
                      m, ...                    % metafile struct, m
                      0, ...                    % 1: if you want to do secondary template matching
                      []);              % filename to store output, leave as [] if you don't want to save
                                                           
%% Extract PD data
m.fps = 180; % (projector frame rate)*3  (*3 for RGB)
% m.pdthr = 3e3; % Can set now, or comment out to set graphically in function (telemetry)
m.pdthr = 3; % OpenEphys projection PD 
[m] = extractTrialADC_PD(filename.raw, ... % Binary File
                        m, ....     % metafile struct, m
                        'test.mat' ); % filename to store output, leave as [] if you don't want to save
                    
%% pd extraction
nLoops = 20;
[m.stimFileName, m.stimFilePath] = uigetfile('*.txt','Select stimGL StimData file:');
stimData=importdata([m.stimFilePath,m.stimFileName]);
m.stimXYPos = stimData.data(:,5:6);

m.StimGL_nloops = nLoops; % Number of loops in stimGL

pdDiff = double(diff(m.pd));
pdDiffThreshold = 1.5e3; % during each repeat of the stimulus there is a repeated PD event, this threshold finds it
[~, m.loopEndIdx]=find(pdDiff>pdDiffThreshold);
m.loopEndIdx = [0, m.loopEndIdx];
if length(m.loopEndIdx) > nLoops
    m.stimLength = diff(m.pd(m.loopEndIdx + 1)); % stimulus length of each repeat
else
    m.stimLength = diff([m.pd(m.loopEndIdx + 1), m.pd(end)]); % stimulus length of each repeat
end
m.loopEndIdx = m.loopEndIdx(1:nLoops);
m.stimLength = m.stimLength(1:nLoops);
stimLengthSubFrames = size(m.stimXYPos,1);

%% raster plot
figure
nLoops = 20;
stimLength = size(m.stimXYPos,1);
spikes_fps=round(s.unit_00*180/30000);  
%end
RasterMask=zeros(1,stimLength);
ind=find(spikes_fps>0 & spikes_fps<stimLength);
RasterMask(spikes_fps(ind))=ones(1,size(ind,2));

T=RasterMask; % linear raster array
for re=2:nLoops
    RasterMask=zeros(1,stimLength);
    temp=spikes_fps-round(m.pd(m.loopEndIdx(re)+1)*180/30000);
    ind=find(temp>0 & temp<stimLength);
    RasterMask(temp(ind))=ones(1,size(ind,2));
    T=[T RasterMask];
end
TE=find(T==1);
subplot(10,1,5:9);  % Spike Raster
    rasterplot(TE,nLoops,stimLength,gca);