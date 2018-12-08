function stim = addProj2M(m, stim, nLoops, pdDiffThreshold)
%% Basic projection info derivation
stim.xMap=stim.xPix/stim.W;
stim.yMap=stim.yPix/stim.H;
dt=1/stim.fps;

%% Get Stimulus Data
[stim.stimFileName, stim.stimFilePath] = uigetfile('*.txt','Select stimGL StimData file:');
stimData=importdata([stim.stimFilePath,stim.stimFileName]);
stim.stimPos = stimData.data(:,5:6);

%% Get phododiode trial start indices and photodiode trial lengths
pdDiff = double(diff(m.pd));
[~, stim.repeatIndex]=find(pdDiff>pdDiffThreshold);

if length(stim.repeatIndex) < nLoops
    stim.repeatIndex(end+1) = m.pd(end);
end
stim.repeatIndex = [0, stim.repeatIndex];
for ii = 1:length(stim.repeatIndex)-1
    stim.stimLength(ii) = m.pd(stim.repeatIndex(ii+1)) - m.pd(stim.repeatIndex(ii)+1);
    stim.trialFrameCounts(ii) = length(find(m.pd >= m.pd(stim.repeatIndex(ii)+1) & m.pd <= m.pd(stim.repeatIndex(ii+1))));
end
stim.repeatIndex = stim.repeatIndex(1:nLoops+1);
stim.stimLength = stim.stimLength(1:nLoops);
stim.trialFrameCounts = stim.trialFrameCounts(1:nLoops);

%% Find number of frames actually shown and find dropped frames
stim.pdIntervals = diff(m.pd(1:stim.repeatIndex(end)));
stim.doubleSkipIdx = m.pd(stim.pdIntervals < pdDiffThreshold & stim.pdIntervals > 3*2.2*stim.sRateHz/stim.fps);
stim.singleSkipIdx = m.pd(stim.pdIntervals <= 3*2.2*stim.sRateHz/stim.fps & stim.pdIntervals > 3*1.2*stim.sRateHz/stim.fps);

%% Reconstruct target pixel and angular data
stim.stimVel=diff(stim.stimPos,1)/dt; % get velocity at each t
stim.stimVel=[stim.stimVel(1,:); stim.stimVel]; % fill in velocity at t = 0 since diff() discards first element

stim.angleStimPos(:,1) = rad2deg(-atan((stim.stimPos(:,1) - stim.xPix/2)/(stim.D*stim.xMap))); % get angular positioning of stimulus
stim.angleStimPos(:,2) = rad2deg(-atan((stim.stimPos(:,2) - stim.yPix/2 - stim.C*stim.yMap)/(stim.D*stim.yMap))); % get angular positioning of stimulus
stim.outOfBoundsIdx = ismember(stim.stimPos,[stim.xPix stim.yPix],'rows') + 1; % get indices of subframes that were not shown by projector
stim.outOfBoundsIdx(find(stim.outOfBoundsIdx == 2)) = nan;
stim.angleStimVel=diff(stim.angleStimPos,1)/dt; % get velocity at each t
stim.angleStimVel=[stim.angleStimVel(1,:); stim.angleStimVel]; % fill in velocity at t = 0 since diff() discards first element

%% Match PD data to stimulus
oobStimPos = stim.stimPos.*stim.outOfBoundsIdx;
oobStimVel = stim.stimVel.*stim.outOfBoundsIdx;
oobAngStimPos = stim.angleStimPos.*stim.outOfBoundsIdx;
oobAngStimVel = stim.angleStimVel.*stim.outOfBoundsIdx;

for ii = 1:length(stim.repeatIndex)-1
    repeatArray = stim.repeatIndex(ii)+1:stim.repeatIndex(ii+1);
    if length(repeatArray) > length(stim.angleStimPos)
        repeatArray = stim.repeatIndex(ii) + 1:stim.repeatIndex(ii) + 1 + length(stim.angleStimPos);
    end
    
    stim.matchPos(repeatArray,:) = stim.stimPos(1:length(repeatArray),:);
    stim.matchVel(repeatArray,:) = stim.stimVel(1:length(repeatArray),:);
    stim.oobMatchPos(repeatArray,:) = oobStimPos(1:length(repeatArray),:);
    stim.oobMatchVel(repeatArray,:) = oobStimVel(1:length(repeatArray),:);
    
    stim.matchAngPos(repeatArray,:) = stim.angleStimPos(1:length(repeatArray),:);
    stim.matchAngVel(repeatArray,:) = stim.angleStimVel(1:length(repeatArray),:);
    stim.oobMatchAngPos(repeatArray,:) = oobAngStimPos(1:length(repeatArray),:);
    stim.oobMatchAngVel(repeatArray,:) = oobAngStimVel(1:length(repeatArray),:);
end
end