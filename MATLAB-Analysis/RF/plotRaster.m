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
stimLength = m.stimLength;
maxStimLength = max(stimLength);
nLoops = m.StimGL_nloops;

if strcmpi(selectUnits, 'all')
    selectUnits = s.clusters;
end

for ii = selectUnits
    if ~isempty(s.(sprintf('unit_%02i',ii)))
        spikeLocations = [];
        for jj=1:nLoops
            %nextTrialShift = s.units{ii} - m.pd(m.repeatIndex(jj) + 1);
            nextTrialShift = double(s.(sprintf('unit_%02i',ii))) - m.pd(m.repeatIndex(jj) + 1);
            nextSpikeTrain = nextTrialShift(nextTrialShift>=0 & nextTrialShift<=stimLength(jj));
            spikeTrain{jj} = nextSpikeTrain;
            spikeLocations = [spikeLocations; nextSpikeTrain];
        end
        
        figure
        title(sprintf('unit_%02i',ii))
        set(gcf,'color','w');
        
        axCellList{1} = subplot(10,1,1:2);  % Velocity angle
        stimAngles = round(atan2d(m.angleStimVel(:,2), m.angleStimVel(:,1)));
        stimAngles(stimAngles == -180) = 180;
        stimAngles(ismember(m.angleStimVel, [0 0], 'rows')) = nan;
        stimAngles(m.outOfBoundsIdx) = nan;
        plot(linspace(0,maxStimLength/m.sRateHz,length(stimAngles-2)), stimAngles,'r');
        xlim([0 maxStimLength/m.sRateHz])
        set(gca,'xtick',[]); set(gca,'ytick',[-180 -90 0 90 180]); grid on
        ylabel('Stimulus trajectory angle (°)'); ylim([-200 200])
        
        axCellList{2} = subplot(10,1,3:8);  % Spike Raster
        rasterplotv2(spikeTrain,maxStimLength,gca,m.sRateHz);
        
        axCellList{3} = subplot(10,1,9:10);   % Spike Histogram
        [nelements,centers]=hist(spikeLocations,0:maxStimLength/80:maxStimLength);
        spikeHis=csaps(centers,nelements,0.5,1:maxStimLength); % smooth and upsample
        plot((1:maxStimLength)/m.sRateHz, spikeHis/nLoops)
        xlim([0 maxStimLength/m.sRateHz])
        xlabel('Time (s)')
        ylabel('Mean spike frequency')
        
        adjustLims(axCellList)
        saveas(gcf, ['raster_unit_' num2str(ii)], 'epsc');
        saveas(gcf, ['raster_unit_' num2str(ii)], 'fig');
    else
        warning(['Unit ' num2str(ii) ' has no spikes. A raster will not be plotted...']);
    end
end

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
