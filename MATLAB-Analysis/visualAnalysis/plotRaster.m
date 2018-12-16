function [] = plotRaster(m, s, stim, selectUnits, opt, saveFig)
%%%% getRaster
%%%% A. Haddad, D. Ko, H. Lin
%%%% INPUT:
%%%% Combine code segments from spikegl_harvest_HTL.m, SpikeRastor.m
%%%% filename.raw=['181017_01.bin']; % Raw data file
%%%% filename.kilosortOutput=['181017_01_sorted.mat']; %including structure s and trial_spikes, trial_cluster

%% Check Required Functions and set up variables
addpath(genpath('requiredFunctions')) % path to folder with functions
reqFuncStr = '';
if ~exist('rasterplotv2.m','file');   reqFuncStr=strcat(reqFuncStr,'rasterplotv2.m, '); end

if ~strcmp(reqFuncStr,'')
    error(['The following functions could not be found in the search path:', newline, reqFuncStr])
end

if strcmpi(selectUnits, 'all')
    selectUnits = 1:length(s.clusters);
end

stimLength = stim.stimLength;
maxStimLength = max(stimLength);
nLoops = stim.StimGL_nloops;
%% spike raster
for ii = selectUnits
    if ~isempty(s.(sprintf('unit_%s',s.clusters(ii))))
        spikeLocations = [];
        
        % get indices of spikes relative to repeatIndex starts and store in
        % a cell array for spike raster and a vector for spike histogram
        for jj=1:nLoops
            %nextTrialShift = s.units{ii} - m.pd(stim.repeatIndex(jj) + 1);
            nextTrialShift = double(s.(sprintf('unit_%s',s.clusters(ii)))) - m.pd(stim.repeatIndex(jj) + 1);
            nextSpikeTrain = nextTrialShift(nextTrialShift>=0 & nextTrialShift<=stimLength(jj));
            spikeTrain{jj} = nextSpikeTrain;
            spikeLocations = [spikeLocations; nextSpikeTrain];
        end
        
        figure
        title(sprintf('unit\\_%s',s.clusters(ii)))
        set(gcf,'color','w');
        
        % plot trajectory angles throughout the experiment, removing points
        % where the trajectory was out of bounds
        axCellList{1} = subplot(10,1,1:2);
        stimAngles = round(atan2d(stim.angleStimVel(:,2), stim.angleStimVel(:,1)));
        stimAngles(stimAngles == -180) = 180;
        stimAngles(ismember(stim.angleStimVel, [0 0], 'rows')) = nan;
        stimAngles(isnan(stim.outOfBoundsIdx)) = nan;
        plot(linspace(0,maxStimLength/m.sRateHz,length(stimAngles-2)), stimAngles,'r');
        xlim([0 maxStimLength/m.sRateHz])
        set(gca,'xtick',[]); set(gca,'ytick',[-180 -90 0 90 180]); grid on
        ylabel('Stimulus trajectory angle (°)'); ylim([-200 200])
        
        % plot spike raster
        axCellList{2} = subplot(10,1,3:8);
        rasterplotv2(spikeTrain,maxStimLength,gca,m.sRateHz*1000);
        
        % plot smoothed and upsampled spike histogram - TO DO: MAKE IT NOT
        % GO NEGATIVE
        axCellList{3} = subplot(10,1,9:10);
        [nelements,centers]=hist(spikeLocations,0:m.sRateHz*opt.rasterBinSize/1000:maxStimLength);
        %spikeHis=csaps(centers,nelements,0.5,1:maxStimLength);
        plot(centers/m.sRateHz, 10*(nelements/nLoops))
        xlim([0 maxStimLength/m.sRateHz])
        xlabel('Time (s)')
        ylabel('Mean spike rate (Hz)')
        
        adjustLims(axCellList)
        
        if saveFig == 2
            export_fig(sprintf('%s_raster_unit_%s.eps',opt.preName,num2str(ii)))
            saveas(gcf, [opt.preName '_raster_unit_' s.clusters(ii)], 'fig');
        elseif saveFig
            saveas(gcf, [opt.preName '_raster_unit_' s.clusters(ii)], 'epsc');
            saveas(gcf, [opt.preName '_raster_unit_' s.clusters(ii)], 'fig');
        end
    else
        warning(['Unit ' s.clusters(ii) ' has no spikes. A raster will not be plotted...']);
    end
end
end

%% UI control
% UI Controls, we need a variable to store the axis handle for each subplot
% we want the axis to change for, additionally each axis should have the
% same x axis scale
function adjustLims(axCellList)
lowerLim = uicontrol(...
    'style',            'edit',...
    'units',            'pixels',...
    'position',         [0 0 80 22],... [LEFT BOTTOM WIDTH HEIGHT]
    'ForegroundColor',  [0 0 0 0],...
    'String',           num2str(axCellList{1}.XLim(1),'%.1e'),...
    'Callback',         {@setAxLim,axCellList,1}...
    );
upperLim = uicontrol(...
    'style',            'edit',...
    'units',            'pixels',...
    'position',         [80 0 80 22],... [LEFT BOTTOM WIDTH HEIGHT]
    'ForegroundColor',  [0 0 0 0],...
    'String',           num2str(axCellList{1}.XLim(2),'%.1e'),...
    'Callback',         {@setAxLim,axCellList,2}...
    );
end

function setAxLim(src, event, axCellList,LorR)
str = get(src, 'String'); % Correct way to get data 'edit' fields
if isnan(str2double(str))     % is not a number?
    set(src, 'String','0');  % reset to 0
    warndlg('Input must be numerical'); % warn dialog
end
xL = axCellList{1}.XLim;
switch LorR
    case 1
        xL(1) = str2double(str);
    case 2
        xL(2) = str2double(str);
end
for ii=1:size(axCellList,2)
    axCellList{ii}.XLim = xL;
end
drawnow
end