function m = addProjectionInfoToM(m, nLoops)
% Basic info
m.D=22.5;        % Distance of the DF head to the screen
m.W=67;       % Width of the projection
m.H=91;        % Height of the projection
m.C=13.5;        % Some distance
m.xPix = 480;
m.yPix = 640;
m.xMap=m.xPix/m.W; % 111; % pix/mm 125mm
m.yMap=m.yPix/m.H;  % 63; % pix/mm % 80mm

% Get Stimulus Data
[m.stimFileName, m.stimFilePath] = uigetfile('*.txt','Select stimGL StimData file:');
stimData=importdata([m.stimFilePath,m.stimFileName]);
m.stimXYPos = stimData.data(:,5:6);

m.StimGL_nloops = nLoops; % Number of loops in stimGL
dt=1/m.fps;    % fps is of projector
pdDiff = double(diff(m.pd));
pdDiffThreshold = 1.5e3; % during each repeat of the stimulus there is a repeated PD event, this threshold finds it
[~, m.repeatIndex]=find(pdDiff>pdDiffThreshold);
m.repeatIndex = [0, m.repeatIndex];
if length(m.repeatIndex) > nLoops
    m.stimLength = diff(m.pd(m.repeatIndex + 1)); % stimulus length of each repeat
else
    m.stimLength = diff([m.pd(m.repeatIndex + 1), m.pd(end)]); % stimulus length of each repeat
end
m.repeatIndex = m.repeatIndex(1:nLoops);
m.stimLength = m.stimLength(1:nLoops);
stimLengthSubFrames = size(m.stimXYPos,1); % ideal stimulus as in stimGL file

%% Reconstruct target angular data
m.angleStimXYPos(:,1) = rad2deg(-atan((m.stimXYPos(:,1) - m.xPix/2)/(m.D*m.xMap))); % get angular positioning of stimulus
m.angleStimXYPos(:,2) = rad2deg(-atan((m.stimXYPos(:,2) - m.yPix/2 - m.C*m.yMap)/(m.D*m.yMap))); % get angular positioning of stimulus
m.outOfBoundsIdx = ismember(m.stimXYPos,[m.xPix m.yPix],'rows');

m.angleStimVel=diff(m.angleStimXYPos,1)/dt; % get velocity at each t
%m.angleStimVel(:,1)=csaps(1:stimLengthSubFrames-1,m.angleStimVel(:,1),0.5,1.5:stimLengthSubFrames-0.5); % smooth out velocity vectors over t
%m.angleStimVel(:,2)=csaps(1:stimLengthSubFrames-1,m.angleStimVel(:,2),0.5,1.5:stimLengthSubFrames-0.5); % smooth out velocity vectors over t
m.angleStimVel=[m.angleStimVel(1,:); m.angleStimVel]; % fill in velocity at t = 0 since diff() discards first element

for ii=1:nLoops
    repeatTime(ii)=round(double(m.pd(m.repeatIndex(ii)+1))*m.fps/m.sRateHz);
    %pixLoopedXYPos = [pixLoopedXYPos; nan((repeatTime(ii)-size(pixLoopedXYPos,1)),2), m.stimXYPos];
    %angleLoopedXYPos = [angleLoopedXYPos; nan((m.repeatTime(ii)*-size(angleLoopedXYPos,1)),2), m.angleStimXYPos];
    %angleLoopedVel = [angleLoopedVel; nan((m.repeatTime(ii)-size(angleLoopedVel,1)),2), m.angleStimVel];
    %m.matchedAngleStimXYPos = 
    %m.matchedAngleStimVel = 
end

end