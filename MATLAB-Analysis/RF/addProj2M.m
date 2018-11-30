function m = addProj2M(m, nLoops, varargin)
nin=nargin;
switch nin
    case 2
        pdDiffThreshold = 1.5e3;
    case 3
        pdDiffThreshold = varargin{1};
    otherwise
        error ('Invalid Arguments');
end

%% Basic projection info
m.D=22.5;               % SET: Distance of the DF head out of the screen
m.W=67;                 % SET: Width of the projection
m.H=91;                 % SET: Height of the projection
m.C=13.5;               % SET: Distance of DF head from centre of screen
m.xPix = 480;           % xaxis pixel count
m.yPix = 640;           % yaxis pixel count

m.xMap=m.xPix/m.W;
m.yMap=m.yPix/m.H;
dt=1/m.fps;

%% Get Stimulus Data
[m.stimFileName, m.stimFilePath] = uigetfile('*.txt','Select stimGL StimData file:');
stimData=importdata([m.stimFilePath,m.stimFileName]);
m.stimXYPos = stimData.data(:,5:6);
m.StimGL_nloops = nLoops; % SET: Number of loops in stimGL

%% Get phododiode trial start indices and photodiode trial lengths
pdDiff = double(diff(m.pd));
[~, m.repeatIndex]=find(pdDiff>pdDiffThreshold);

if length(m.repeatIndex) < nLoops
    m.repeatIndex(end+1) = m.pd(end);
end
m.repeatIndex = [0, m.repeatIndex];
for ii = 1:length(m.repeatIndex)-1
    m.stimLength(ii) = m.pd(m.repeatIndex(ii+1)) - m.pd(m.repeatIndex(ii)+1);
    m.trialFrameCounts(ii) = length(find(m.pd >= m.pd(m.repeatIndex(ii)+1) & m.pd <= m.pd(m.repeatIndex(ii+1))));
end
m.repeatIndex = m.repeatIndex(1:nLoops+1);
m.stimLength = m.stimLength(1:nLoops);
m.trialFrameCounts = m.trialFrameCounts(1:nLoops);

%% Find number of frames actually shown and find dropped frames
m.pdIntervals = diff(m.pd(1:m.repeatIndex(end)));
m.doubleSkipIdx = m.pd(m.pdIntervals < pdDiffThreshold & m.pdIntervals > 3*2.2*m.sRateHz/m.fps);
m.singleSkipIdx = m.pd(m.pdIntervals <= 3*2.2*m.sRateHz/m.fps & m.pdIntervals > 3*1.2*m.sRateHz/m.fps);

%% Reconstruct target angular data
m.angleStimXYPos(:,1) = rad2deg(-atan((m.stimXYPos(:,1) - m.xPix/2)/(m.D*m.xMap))); % get angular positioning of stimulus
m.angleStimXYPos(:,2) = rad2deg(-atan((m.stimXYPos(:,2) - m.yPix/2 - m.C*m.yMap)/(m.D*m.yMap))); % get angular positioning of stimulus
m.outOfBoundsIdx = ismember(m.stimXYPos,[m.xPix m.yPix],'rows'); % get indices of subframes that were not shown by projector

m.angleStimVel=diff(m.angleStimXYPos,1)/dt; % get velocity at each t
m.angleStimVel=[m.angleStimVel(1,:); m.angleStimVel]; % fill in velocity at t = 0 since diff() discards first element

%% Match PD data to stimulus
for ii = 1:length(m.repeatIndex)-1
    repeatArray = m.repeatIndex(ii)+1:m.repeatIndex(ii+1);
    if length(repeatArray) > length(m.angleStimXYPos)
        repeatArray = m.repeatIndex(ii) + 1:m.repeatIndex(ii) + 1 + length(m.angleStimXYPos);
    end
    m.matchedAngleStimXYPos(repeatArray,:) = m.angleStimXYPos(1:length(repeatArray),:);
    m.matchedAngleStimVel(repeatArray,:) = m.angleStimVel(1:length(repeatArray),:);
end
end