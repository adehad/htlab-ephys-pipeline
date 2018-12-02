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
m.stimPos = stimData.data(:,5:6);
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

%% Reconstruct target pixel and angular data
m.stimVel=diff(m.stimPos,1)/dt; % get velocity at each t
m.stimVel=[m.stimVel(1,:); m.stimVel]; % fill in velocity at t = 0 since diff() discards first element

m.angleStimPos(:,1) = rad2deg(-atan((m.stimPos(:,1) - m.xPix/2)/(m.D*m.xMap))); % get angular positioning of stimulus
m.angleStimPos(:,2) = rad2deg(-atan((m.stimPos(:,2) - m.yPix/2 - m.C*m.yMap)/(m.D*m.yMap))); % get angular positioning of stimulus
m.outOfBoundsIdx = ismember(m.stimPos,[m.xPix m.yPix],'rows') + 1; % get indices of subframes that were not shown by projector
m.outOfBoundsIdx(find(m.outOfBoundsIdx == 2)) = nan;
m.angleStimVel=diff(m.angleStimPos,1)/dt; % get velocity at each t
m.angleStimVel=[m.angleStimVel(1,:); m.angleStimVel]; % fill in velocity at t = 0 since diff() discards first element

%% Match PD data to stimulus
oobStimPos = m.stimPos.*m.outOfBoundsIdx;
oobStimVel = m.stimVel.*m.outOfBoundsIdx;
oobAngStimPos = m.angleStimPos.*m.outOfBoundsIdx;
oobAngStimVel = m.angleStimVel.*m.outOfBoundsIdx;

for ii = 1:length(m.repeatIndex)-1
    repeatArray = m.repeatIndex(ii)+1:m.repeatIndex(ii+1);
    if length(repeatArray) > length(m.angleStimPos)
        repeatArray = m.repeatIndex(ii) + 1:m.repeatIndex(ii) + 1 + length(m.angleStimPos);
    end
    
    m.matchPos(repeatArray,:) = m.stimPos(1:length(repeatArray),:);
    m.matchVel(repeatArray,:) = m.stimVel(1:length(repeatArray),:);
    m.oobMatchPos(repeatArray,:) = oobStimPos(1:length(repeatArray),:);
    m.oobMatchVel(repeatArray,:) = oobStimVel(1:length(repeatArray),:);
    
    m.matchAngPos(repeatArray,:) = m.angleStimPos(1:length(repeatArray),:);
    m.matchAngVel(repeatArray,:) = m.angleStimVel(1:length(repeatArray),:);
    m.oobMatchAngPos(repeatArray,:) = oobAngStimPos(1:length(repeatArray),:);
    m.oobMatchAngVel(repeatArray,:) = oobAngStimVel(1:length(repeatArray),:);
end
end