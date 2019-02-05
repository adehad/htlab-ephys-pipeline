function stim = addProj2Stim(m, stim, pdDiffThreshold, pType)
%% Basic projection info derivation
if strcmpi(pType, 'flat')
    stim.xMap=stim.xPix/stim.W;
    stim.yMap=stim.yPix/stim.H;
elseif strcmpi(pType, 'dome')
    
end
dt=1/m.fps;

%% Get stimulus data
[stim.stimFileName, stim.stimFilePath] = uigetfile('*.txt','Select stimGL StimData file:');
stimData=importdata([stim.stimFilePath,stim.stimFileName]);
stim.stimPos = stimData.data(:,5:6);

%% Get photodiode data

% find loop start and end indices in m.pd
pdDiff = double(diff(m.pd));
[~, stim.loopEndIdx]=find(pdDiff>pdDiffThreshold);

if length(stim.loopEndIdx) < stim.nLoopsPre
    stim.loopEndIdx(end+1) = length(m.pd);
    disp('Using final pd event as endtime of last loop...')
elseif length(stim.loopEndIdx) > stim.nLoopsPre
    stim.loopEndIdx = stim.loopEndIdx(1:stim.nLoopsPre); % cut post-experiment artifacts
    disp('Removing any detected loopEndIdx after nLoopsPre...')
end

if length(stim.loopEndIdx) == stim.nLoopsPre
    disp('Correct number of loopEndIdx detected!')
else
    error('Number of loopEndIdx does not equal nLoopsPre!')
end

stim.loopStartIdx = [1 stim.loopEndIdx(1:end-1) + 1];

% find stimulus lengths of each loop in time and in frames
for ii = 1:stim.nLoopsPre
    stim.stimLength(ii) = m.pd(stim.loopEndIdx(ii)) - m.pd(stim.loopStartIdx(ii));
    stim.frameCount(ii) = stim.loopEndIdx(ii)- stim.loopStartIdx(ii);
end

%% Discard bad loops
goodL = (stim.stimLength < 1.01*mode(stim.stimLength)) & (stim.stimLength > 0.99*mode(stim.stimLength));
stim.nLoopsPost = nnz(goodL);
if stim.nLoopsPost ~= stim.nLoopsPre
    stim.loopStartIdx = stim.loopStartIdx(goodL);
    stim.loopEndIdx = stim.loopEndIdx(goodL);
    stim.stimLength = stim.stimLength(goodL);
    stim.frameCount = stim.frameCount(goodL);
    warning('Loops [%s] are too long/short in time... discarding...', dispArray(find(goodL)));
end

%% Reconstruct target pixel data and find out of bounds indices
if strcmpi(pType, 'flat')
    stim.stimVel=diff(stim.stimPos,1)/dt; % get velocity at each t
    stim.stimVel=[stim.stimVel(1,:); stim.stimVel]; % fill in velocity at t = 0 since diff() discards first element

    stim.outOfBoundsIdx = ismember(stim.stimPos,[stim.xPix stim.yPix],'rows'); % get indices of subframes that were not shown by projector
    stim.outOfBoundsIdx(stim.outOfBoundsIdx == 1) = nan;
end

%% Reconstruct target angle data
if strcmpi(pType, 'flat')
    stim.angleStimPos(:,1) = rad2deg(-atan((stim.stimPos(:,1) - stim.xPix/2 - stim.C*stim.yMap)/(stim.E*stim.xMap))); % get angular positioning of stimulus
    stim.angleStimPos(:,2) = rad2deg(-atan((stim.stimPos(:,2) - stim.yPix/2 - stim.D*stim.yMap)/(stim.E*stim.yMap))); % get angular positioning of stimulus

    stim.angleStimVel=diff(stim.angleStimPos,1)/dt; % get velocity at each t
    stim.angleStimVel=[stim.angleStimVel(1,:); stim.angleStimVel]; % fill in velocity at t = 0 since diff() discards first element
    anom = union(find(abs(stim.angleStimVel(:,1))>150),find(abs(stim.angleStimVel(:,2))>150));
    stim.angleStimVel(anom,:) = nan;
end

%% Match PD data to stimulus
for ii = 1:length(stim.loopEndIdx)-1
    temp = stim.loopStartIdx(ii):stim.loopEndIdx(ii);
    if length(temp) > length(stim.stimPos)/3
        temp = temp(1:(length(stim.stimPos)/3));
    end
    dSamp = 1:3:3*length(temp);
    
    stim.matchPos(temp,:) = stim.stimPos(dSamp,:);
    stim.matchVel(temp,:) = stim.stimVel(dSamp,:);
    
    stim.matchAngPos(temp,:) = stim.angleStimPos(dSamp,:);
    stim.matchAngVel(temp,:) = stim.angleStimVel(dSamp,:);
end

end