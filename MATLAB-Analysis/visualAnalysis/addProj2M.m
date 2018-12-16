function stim = addProj2M(m, stim, pdDiffThreshold, correctSpikes)
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
stim.repeatIndex = [0, stim.repeatIndex]; %append 0 as end of zeroth trial
stim.repeatIndex = stim.repeatIndex(1:stim.StimGL_nloops+1);
for ii = 1:length(stim.repeatIndex)-1
    stim.stimLength(ii) = m.pd(stim.repeatIndex(ii+1)) - m.pd(stim.repeatIndex(ii)+1);
    stim.trialFrameCounts(ii) = length(find(m.pd >= m.pd(stim.repeatIndex(ii)+1) & m.pd <= m.pd(stim.repeatIndex(ii+1))));
end
% stim.repeatIndex = stim.repeatIndex(1:stim.StimGL_nloops+1);
% stim.stimLength = stim.stimLength(1:stim.StimGL_nloops);
% stim.trialFrameCounts = stim.trialFrameCounts(1:stim.StimGL_nloops);

%% Find number of frames actually shown and find obvious dropped frames
stim.pdIntervals = diff(m.pd(1:stim.repeatIndex(end)));
stim.doubleSkipIdx = m.pd(stim.pdIntervals < pdDiffThreshold & stim.pdIntervals > 3*2.2*m.sRateHz/m.fps);
stim.singleSkipIdx = m.pd(stim.pdIntervals <= 3*2.2*m.sRateHz/m.fps & stim.pdIntervals > 3*1.2*m.sRateHz/m.fps);

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

oobAngStimPos = stim.angleStimPos.*stim.outOfBoundsIdx;
oobAngStimVel = stim.angleStimVel.*stim.outOfBoundsIdx;

%% Match PD data to stimulus

for ii = 1:length(stim.repeatIndex)-1
    temp = m.pd(stim.repeatIndex(ii)+1):m.pd(stim.repeatIndex(ii+1));
    if m.trialFrameCounts(ii) > length(stim.angleStimPos)
        temp = stim.repeatIndex(ii) + 1:stim.repeatIndex(ii) + 1 + length(stim.angleStimPos);
    elseif m.trialFrameCounts(ii) < length(stim.angleStimPos)
        
    end
    
    stim.matchPos(temp,:) = stim.stimPos(1:length(temp),:);
    stim.matchVel(temp,:) = stim.stimVel(1:length(temp),:);
    stim.oobMatchPos(temp,:) = oobStimPos(1:length(temp),:);
    stim.oobMatchVel(temp,:) = oobStimVel(1:length(temp),:);
    
    stim.matchAngPos(temp,:) = stim.angleStimPos(1:length(temp),:);
    stim.matchAngVel(temp,:) = stim.angleStimVel(1:length(temp),:);
    stim.oobMatchAngPos(temp,:) = oobAngStimPos(1:length(temp),:);
    stim.oobMatchAngVel(temp,:) = oobAngStimVel(1:length(temp),:);
end

for ii = 1:length(stim.repeatIndex)-1
    temp = stim.repeatIndex(ii)+1:stim.repeatIndex(ii+1);
    if length(temp) > length(stim.angleStimPos)
        temp = stim.repeatIndex(ii) + 1:stim.repeatIndex(ii) + 1 + length(stim.angleStimPos);
    elseif length(temp) < length(stim.angleStimPos)
        
    end
    
    stim.matchPos(temp,:) = stim.stimPos(1:length(temp),:);
    stim.matchVel(temp,:) = stim.stimVel(1:length(temp),:);
    stim.oobMatchPos(temp,:) = oobStimPos(1:length(temp),:);
    stim.oobMatchVel(temp,:) = oobStimVel(1:length(temp),:);
    
    stim.matchAngPos(temp,:) = stim.angleStimPos(1:length(temp),:);
    stim.matchAngVel(temp,:) = stim.angleStimVel(1:length(temp),:);
    stim.oobMatchAngPos(temp,:) = oobAngStimPos(1:length(temp),:);
    stim.oobMatchAngVel(temp,:) = oobAngStimVel(1:length(temp),:);
end

disp('Done loading stimulus data!')

in each repeatIndexSlice

find m.pd slice time length
find true m.pd time length
find difference in time length as percentage
times spike times by percentage
output 
end