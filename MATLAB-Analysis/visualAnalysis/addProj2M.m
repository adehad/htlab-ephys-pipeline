function stim = addProj2M(m, stim, pdDiffThreshold)
%% Basic projection info derivation
stim.xMap=stim.xPix/stim.W;
stim.yMap=stim.yPix/stim.H;
dt=1/m.fps;

%% Get Stimulus Data
[stim.stimFileName, stim.stimFilePath] = uigetfile('*.txt','Select stimGL StimData file:');
stimData=importdata([stim.stimFilePath,stim.stimFileName]);
stim.stimPos = stimData.data(:,5:6);

%% Get phododiode trial end indices and photodiode trial lengths
pdDiff = double(diff(m.pd));
[~, stim.repeatIndex]=find(pdDiff>pdDiffThreshold);

if length(stim.repeatIndex) < stim.StimGL_nloops
    stim.repeatIndex(end+1) = length(m.pd);
end
stim.repeatIndex = stim.repeatIndex(1:stim.StimGL_nloops); % cut post-experiment artifacts

stim.repeatIndex = [0 stim.repeatIndex];

for ii = 1:length(stim.repeatIndex)-1
    stim.stimLength(ii) = m.pd(stim.repeatIndex(ii+1)) - m.pd(stim.repeatIndex(ii)+1);
    stim.frameCount(ii) = stim.repeatIndex(ii+1)- stim.repeatIndex(ii)+1;
end
% stim.repeatIndex = stim.repeatIndex(1:stim.StimGL_nloops+1);
% stim.stimLength = stim.stimLength(1:stim.StimGL_nloops);
% stim.frameCount = stim.frameCount(1:stim.StimGL_nloops);

%% Discard bad trials
goodTrials = (stim.stimLength < 1.01*mode(stim.stimLength)) & (stim.stimLength > 0.99*mode(stim.stimLength));
stim.repeatIndex = [stim.repeatIndex(goodTrials) stim.repeatIndex(end)];
stim.stimLength = stim.stimLength(goodTrials);
stim.frameCount = stim.frameCount(goodTrials);

%% Find number of frames actually shown and find obvious dropped frames
% stim.pdIntervals = diff(m.pd(1:stim.repeatIndex(end)));
% stim.doubleSkipIdx = m.pd(stim.pdIntervals < pdDiffThreshold & stim.pdIntervals > 3*2.2*m.sRateHz/m.fps);
% stim.singleSkipIdx = m.pd(stim.pdIntervals <= 3*2.2*m.sRateHz/m.fps & stim.pdIntervals > 3*1.2*m.sRateHz/m.fps);

%% Reconstruct target pixel data and find out of bounds indices
stim.stimVel=diff(stim.stimPos,1)/dt; % get velocity at each t
stim.stimVel=[stim.stimVel(1,:); stim.stimVel]; % fill in velocity at t = 0 since diff() discards first element

stim.outOfBoundsIdx = ismember(stim.stimPos,[stim.xPix stim.yPix],'rows') + 1; % get indices of subframes that were not shown by projector
stim.outOfBoundsIdx(find(stim.outOfBoundsIdx == 2)) = nan;

oobStimPos = stim.stimPos.*stim.outOfBoundsIdx;
oobStimVel = stim.stimVel.*stim.outOfBoundsIdx;

%% Reconstruct target angle data
stim.angleStimPos(:,1) = rad2deg(-atan((stim.stimPos(:,1) - stim.xPix/2 - stim.C*stim.yMap)/(stim.E*stim.xMap))); % get angular positioning of stimulus
stim.angleStimPos(:,2) = rad2deg(-atan((stim.stimPos(:,2) - stim.yPix/2 - stim.D*stim.yMap)/(stim.E*stim.yMap))); % get angular positioning of stimulus


stim.angleStimVel=diff(stim.angleStimPos,1)/dt; % get velocity at each t
stim.angleStimVel=[stim.angleStimVel(1,:); stim.angleStimVel]; % fill in velocity at t = 0 since diff() discards first element
anom = union(find(abs(stim.angleStimVel(:,1))>150),find(abs(stim.angleStimVel(:,2))>150));
stim.angleStimVel(anom,:) = nan;

% discont = abs(diff(stim.angleStimPos)) > 8;
% trajEnds = find(discont(:,1) | discont(:,2));
% trajEnds(1:2:end) = trajEnds(1:2:end)+1;
% stim.angleStimVel = zeros(size(stim.angleStimPos));
% for ii = 1:2:length(trajEnds)
%     trTime = trajEnds(ii+1)-trajEnds(ii)+1;
%     sp = m.fps*(stim.angleStimPos(trajEnds(ii+1),:)-stim.angleStimPos(trajEnds(ii),:))/trTime;
%     stim.angleStimVel(trajEnds(ii):trajEnds(ii+1),:) = repmat(sp,trajEnds(ii+1)-trajEnds(ii)+1,1);
%     stim.trajMatch(trajEnds(ii):trajEnds(ii+1),1) = ii;
%     stim.trajMatch(trajEnds(ii):trajEnds(ii+1),2) = (1:length(trajEnds(ii):trajEnds(ii+1)))/length(trajEnds(ii):trajEnds(ii+1));
% end

% stim.trajMatch = cat(1,stim.trajMatch,zeros(length(stim.angleStimVel)-length(stim.trajMatch),2));

oobAngStimPos = stim.angleStimPos.*stim.outOfBoundsIdx;
oobAngStimVel = stim.angleStimVel.*stim.outOfBoundsIdx;

%% Match PD data to stimulus

for ii = 1:length(stim.repeatIndex)-1
    temp = stim.repeatIndex(ii)+1:stim.repeatIndex(ii+1);
    if length(temp) > length(stim.stimPos)/3
        temp = temp(1:(length(stim.stimPos)/3));
    end
    dSamp = 1:3:3*length(temp);
    
    %stim.pdTraj(temp,:) = stim.trajMatch(dSamp,:);
    stim.matchPos(temp,:) = stim.stimPos(dSamp,:);
    stim.matchVel(temp,:) = stim.stimVel(dSamp,:);
    stim.oobMatchPos(temp,:) = oobStimPos(dSamp,:);
    stim.oobMatchVel(temp,:) = oobStimVel(dSamp,:);
    
    stim.matchAngPos(temp,:) = stim.angleStimPos(dSamp,:);
    stim.matchAngVel(temp,:) = stim.angleStimVel(dSamp,:);
    stim.oobMatchAngPos(temp,:) = oobAngStimPos(dSamp,:);
    stim.oobMatchAngVel(temp,:) = oobAngStimVel(dSamp,:);
end

disp(['Done loading stimulus data! Trials discarded: [1 ' num2str(find(~goodTrials)+1) ']'])
end