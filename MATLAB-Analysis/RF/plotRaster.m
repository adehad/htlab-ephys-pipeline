function [] = getRaster(m, s, selectUnits)
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

%% set up some parameters in structure m
% [m, fpath, mfile] = readMetafile;
% m.metafile = mfile;
% m.metapath = fpath;
% NO MetaFile with Daniel's Data, manually entering:
m.StimGL_nloops = 20;

%% start to construct data for rastor plot
% Get Stimulus Data
[m.stimFileName, m.stimFilePath] = uigetfile('*.txt','Select stimGL StimData file:');
stimData=importdata([m.stimFilePath,m.stimFileName]);
m.stimXYPos = stimData.data(:,5:6); %int16(xytraj);

dt=1/m.fps;    % 1/fps    fps is of projector
stimLength=size(m.stimXYPos,1);

for ii = length(s.units)
    alignedUnits{ii} = double(s.units{ii} - m.pd(1) + 1); %Start aligned the spikes  /m.msec;
end
pdDiff = double(diff(m.pd));

%% Step 1: Reconstruct target angular data
m.angleStimXYPos(:,1) = rad2deg(-atan((m.stimXYPos(:,1) - xPix/2)/(D*xMap)));
m.angleStimXYPos(:,2) = rad2deg(-atan((m.stimXYPos(:,2) - yPix/2 - C*yMap)/(D*yMap)));

m.angleStimXYVel=diff(m.angleStimXYPos,1)/dt;
m.angleStimXYVel(:,1)=csaps(1:stimLength-1,m.angleStimXYVel(:,1),0.5,1.5:stimLength-0.5);%?
m.angleStimXYVel(:,2)=csaps(1:stimLength-1,m.angleStimXYVel(:,2),0.5,1.5:stimLength-0.5);%?
m.angleStimXYVel=[m.angleStimXYVel(1,:); m.angleStimXYVel];

%% Step 2: Align all the data series
pdDiffThreshold = 1.6e3; % during each repeat of the stimulus there is a repeated PD event, this threshold finds it
[~, repeatIndex]=find(pdDiff>pdDiffThreshold);
nLoops = m.StimGL_nloops;
pixTargetTraj = m.stimXYPos;
% StartFs =StartFrame;
angleTargetTraj=m.angleStimXYPos;
angleTargetVel=m.angleStimXYVel;
for ii=1:nLoops-1
    repeatTime(ii)=round(double(m.pd(repeatIndex(ii)+1)-m.pd(1))*m.fps/m.sRateHz);
    pixTargetTraj = [pixTargetTraj; zeros((repeatTime(ii)-size(pixTargetTraj,1)),2)];
    pixTargetTraj = [pixTargetTraj; m.stimXYPos];
%     StartFs = [StartFs RepeatTime(nl).*ones(1,size(StartFrame,1))+StartFrame];    
    angleTargetTraj = [angleTargetTraj; zeros((repeatTime(ii)*-size(angleTargetTraj,1)),2)];
    angleTargetTraj = [angleTargetTraj; m.angleStimXYPos];
    angleTargetVel = [angleTargetVel; zeros((repeatTime(ii)-size(angleTargetVel,1)),2)];
    angleTargetVel = [angleTargetVel; m.angleStimXYVel];    
end

%% Step 3: Compile spike raster

if strcmpi(selectUnits, 'all')
    selectUnits = length(alignedUnits);
end
for ii = selectUnits
    singleUnit=round(alignedUnits{ii}*m.fps/m.sRateHz);
    rasterMask=zeros(1,stimLength);
    ind=find(singleUnit>0 & singleUnit<stimLength);
    rasterMask(singleUnit(ind))=ones(1,size(ind,2));

raster(1,:)=rasterMask; 
spikeTiming=find(raster(1,:));
T=rasterMask; % linear raster array
for jj=1:nLoops-1
    rasterMask=zeros(1,stimLength);
    temp=singleUnit-repeatTime(jj);
    ind=find(temp>0 & temp<stimLength);
    rasterMask(temp(ind))=ones(1,size(ind,2));
%     SpikeTiming=[SpikeTiming temp(ind)];
    spikeTiming=[spikeTiming find(rasterMask)];
    raster(jj+1,:)=rasterMask; 
    T=[T rasterMask];
end
TE=find(T==1);
% SpikeTiming=spikes_360(ind);
[nelements,centers]=hist(spikeTiming,0:36:stimLength);
spikeHis=csaps(centers,nelements,0.5,1:stimLength); % smooth and upsample

figure
subplot(10,1,1:2); % Position
    plot(m.angleStimXYPos(:,1),'r'); hold on; plot(m.angleStimXYPos(:,2),'b');
    xlim([0 stimLength])
    ylim([-60 60])
subplot(10,1,3:4);  % Velocity
    plot(m.angleStimXYVel(:,1),'r'); hold on; plot(m.angleStimXYVel(:,2),'b');
    xlim([0 stimLength])
    ylim([-150 150])

if nLoops>1
    ntrials = nLoops;
else
    ntrials = nLoops;
end
subplot(10,1,5:9);  % Spike Raster
    rasterplot(TE,ntrials,stimLength,gca);
%     ylabel(['Trails (',num2str(ntrials),')']) % incorporated to rasterplot
    % rasterplot(TE,nloops-1,stimLength,'plotwidth',0.5);
subplot(10,1,10);   % Spike Hist
    %bar(centers./360,nelements);
    plot((1:stimLength)./360,spikeHis)
    xlim([0 stimLength./360])
    ylim([0 50])
end
%% Plot Receptive Fields
%%%% INPUT: TarTraj (above), Gaze_rest, head_offset, spikes_180 (spike times in 180Hz projector indices)
%%%% Requires: plotellipse.m
%% set up some variables
head_offset_ang=0;  %this depends on the head marker mounting azimuthal accuracy per animal
Gaze_rest=51;       %resting gaze elevation in degree
spikes_180 = singleUnit';
rejected_spikes_180 = spikes_180(spikes_180 <= 0 );
spikes_180 = spikes_180(spikes_180 > 0 );

disp(['# Negative spikes omitted = ', num2str(length(rejected_spikes_180))]);


% saccade cone parameters:
ab=44*pi/180; %88
bb=56.5*pi/180; %113/2=
zb=[0 75*pi/180];%estimated saccade cone axis 75deg
alphab=0;

% takeoff cone parameters:
ab2=22.5*pi/180; %45
bb2=28.5*pi/180; %57
zb2=[0 80*pi/180];%estimated saccade cone axis 80deg
alphab2=0;
spike_tail = 5; % how many samples to plot leading up to the spike location
response_offset = 5; % how many samples to shift spike location in receptive field to account for animal response time
%% plot spike time stimuli in angular coordinate
figure(3); hold on;
ax_fig3 = gca;
plot(m.angleStimXYPos(:,1), m.angleStimXYPos(:,2)+Gaze_rest+head_offset_ang,'.','color',[0.8 0.8 0.8]); 
errorInLoop = [];

% Fix error in KiloSort - can comment out once fixed
spikes_180 = unique(spikes_180);

for ss=1:ntrials
    
    spikes_inTrajTrial = ( spikes_180>( (ss-1)*(size(angleStimXYPos,1)) ) ) & ( spikes_180<( ss*(size(angleStimXYPos,1)) ) ) ;
       
    spikes_inTrajTrial = spikes_180(spikes_inTrajTrial); %
    
    spikes_inTrajTrial = spikes_inTrajTrial - response_offset;
    
    spikes_inTraj.(['trial_',num2str(ss)]) = spikes_inTrajTrial ;
        % Stores what spikes were found for a given repeat of stimuli
        
    
    fprintf('Plotting Trial: %02d, %03d (Spikes) \n', ss, length(spikes_inTrajTrial));
    title(['Plotting Trials 1 to ' , num2str(ss)])
    for kk=1:length(spikes_inTrajTrial)
        temp_startElem = spikes_inTrajTrial(kk);
        temp_tailElem = temp_startElem - spike_tail;
        if temp_tailElem < (ss-1)*(size(angleStimXYPos,1))
            temp_tailElem = (ss-1)*(size(angleStimXYPos,1))+1;
        end
        
%--------% Spike Location
         plot(angleTargetTraj(temp_startElem,1), ...
              angleTargetTraj(temp_startElem,2)+Gaze_rest+head_offset_ang,'r.');
     
%--------% Spike tail
        xTrajVec = angleTargetTraj( temp_tailElem:temp_startElem, 1 );
        yTrajVec = angleTargetTraj( temp_tailElem:temp_startElem, 2 )+Gaze_rest+head_offset_ang;
        
        % filter out any movement that is too large
            movementThresh = 12; % in angular degrees
            
            validVectorElems = abs(diff(xTrajVec))<movementThresh & abs(diff(yTrajVec))<movementThresh;

            temp_EndElem = find(validVectorElems==true,1,'last');
            if isempty(temp_EndElem) % empty if there are no valid elements
                temp_EndElem = length(validVectorElems);
            end

            temp_StartElem = find(validVectorElems(1:temp_EndElem)==false,1,'last');
            if isempty(temp_StartElem) % empty if all remaining are valid
                temp_StartElem = 0; % will +1 below
            end

            validVectorElems = (temp_StartElem+1:temp_EndElem+1);
        try
        xTrajVec = xTrajVec(validVectorElems);
        yTrajVec = yTrajVec(validVectorElems);
        catch ME
            warning(ME.message)
        end
%         plot( ax_fig3, xTrajVec, yTrajVec, 'b' ); % Faster but less useful I think
        plot_dir_ade( xTrajVec, ...            
                      yTrajVec, [], ax_fig3);
%           drawnow;   % use the one outside this for loop to be faster  
    end
    drawnow;  % use this if you want to see the spikes per trial live
end

 
plotellipse(zb*180/pi, bb*180/pi, ab*180/pi, alphab, 'c--')
plot(zb(1)*180/pi,zb(2)*180/pi,'c+')

plotellipse(zb2*180/pi, bb2*180/pi, ab2*180/pi, alphab2, 'm--')
plot(zb2(1)*180/pi,zb2(2)*180/pi,'m+')

plot(0,Gaze_rest,'k+'); axis equal
    xlim([-60 60]); ylim([0 120])
    xlabel('azimuth (deg)'); ylabel('elevation (deg)'); 
    title('Target in body oriented global ref')
figure(3)
hold off
%% plot spike time stimuli in angular coordinate [Heat Map]
% Using TarTraj_s - as this includes the repeats
figure(4); hold on;
% scatter(TarTraj(:,1), TarTraj(:,2)+Gaze_rest+head_offset_ang,10,[0.8 0.8 0.8],'filled'); 
% alpha 0.5
numBins = 6; % Number of bins for each dimensions
histogram2( angleTargetTraj(spikes_180(spikes_180<=length(angleTargetTraj)),1), ...
            angleTargetTraj(spikes_180(spikes_180<=length(angleTargetTraj)),2)+Gaze_rest+head_offset_ang, ...
            numBins, ...
            'DisplayStyle','tile','ShowEmptyBins','on' )
% Omits spikes that are outside the length of the stimuli

plotellipse(zb*180/pi, bb*180/pi, ab*180/pi, alphab, 'c--')
plot(zb(1)*180/pi,zb(2)*180/pi,'c+')

plotellipse(zb2*180/pi, bb2*180/pi, ab2*180/pi, alphab2, 'm--')
plot(zb2(1)*180/pi,zb2(2)*180/pi,'m+')

plot(0,Gaze_rest,'k+'); axis equal
    xlim([-60 60]); ylim([0 120])
    xlabel('azimuth (deg)'); ylabel('elevation (deg)'); 
    title('Target in body oriented global ref')

%% overlay a head map of spike density
