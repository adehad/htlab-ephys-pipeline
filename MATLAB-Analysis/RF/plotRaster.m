function [] = plotRaster(m, s, selectUnits)
%%%% getRaster
%%%% A. Haddad, D. Ko, H. Lin
%%%% INPUT:
%%%% Combine code segments from spikegl_harvest_HTL.m, SpikeRastor.m
%%%% filename.raw=['181017_01.bin']; % Raw data file
%%%% filename.kilosortOutput=['181017_01_sorted.mat']; %including structure s and trial_spikes, trial_cluster

%% Check Required Functions
addpath(genpath('requiredFunctions')) % path to folder with functions
reqFuncStr = '';
if ~exist('readMetafile.m','file'); reqFuncStr=strcat(reqFuncStr,'readMetafile.m, ');           end
if ~exist('get_ch_thr.m','file');   reqFuncStr=strcat(reqFuncStr,'get_ch_thr.m, ');             end
if ~exist('get_events.m','file');   reqFuncStr=strcat(reqFuncStr,'get_events.m (get_ch_thr), ');end
if ~exist('splitconv.m','file');    reqFuncStr=strcat(reqFuncStr,'splitconv.m (get_events), '); end
if ~exist('tconv.m','file');        reqFuncStr=strcat(reqFuncStr,'tconv.m (splitconv), ');      end
if ~exist('wavecull.m','file');     reqFuncStr=strcat(reqFuncStr,'wavecull.m, ');               end
if ~exist('rasterplot.m','file');   reqFuncStr=strcat(reqFuncStr,'rasterplot.m, ');             end
if ~exist('plot_dir_ade.m','file');  reqFuncStr=strcat(reqFuncStr,'plot_dir_ade.m (can replace with plot), ');            end
if ~exist('plotellipse.m','file');  reqFuncStr=strcat(reqFuncStr,'plotellipse.m, ');            end

if ~strcmp(reqFuncStr,'')
    error(['The following functions could not be found in the search path:', newline, reqFuncStr])
end

%% Step 3: Compile spike raster
if strcmpi(selectUnits, 'all')
    selectUnits = length(s.units);
end
for ii = selectUnits
    spikeTrain = [];
    for jj=1:nLoops
        nextTrialShift = s.units{ii} - m.pd(m.repeatIndex + 1);
        nextSpikeTrain = nextTrialShift(find(nextTrialShift>0 & nextTrialShift<m.stimLength));
        spikeTrain=[spikeTrain nextSpikeTrain];
    end
    
    figure
    hold on;
    subplot(10,1,1:2); % Position
    plot(m.angleStimXYPos(:,1),'r');
    plot(m.angleStimXYPos(:,2),'b');
    xlim([0 m.stimLength])
    ylim([-max(m.xPix, m.yPix)/2, max(m.xPix, m.yPix)/2])

%     subplot(10,1,3:4);  % Velocity
%     plot(m.angleStimXYVel(:,1),'r'); hold on; plot(m.angleStimXYVel(:,2),'b');
%     xlim([0 m.stimLength])
%     ylim([-150 150])

    subplot(10,1,5:9);  % Spike Raster
    rasterplot(spikeTrain,nLoops,m.stimLength,gca);
    
    subplot(10,1,10);   % Spike Histogram
    [nelements,centers]=hist(spikeTrain,0:36:m.stimLength);
    spikeHis=csaps(centers,nelements,0.5,1:m.stimLength); % smooth and upsample
    plot((1:m.stimLength)./m.fps,spikeHis)
    xlim([0 m.stimLength./m.fps])
    ylim([0 50])
end