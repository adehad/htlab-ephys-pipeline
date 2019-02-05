function lastFrameIdxList = spikeTime2FrameTime(m, s, stim, selectUnits, opt)
%SPIKETIME2FRAMETIME Summary of this function goes here
%   Detailed explanation goes here
for ii = selectUnits
    singleUnit = double(s.(sprintf('unit_%s',s.clusters(ii)))); % get sing unit spikes
    clear sqC sqMean
    if ~isempty(singleUnit)
        % get spikes for full experiment
        singleUnit = singleUnit(m.pd(1) <= singleUnit & singleUnit <= m.pd(stim.loopEndIdx(end)));
        
        % shift m.pd by the times of each spike to find out the index of last pd
        % event before the spike occurs. Introduce dragonfly neural latency
        % if wanted
        lastFrameIdx = ones(size(singleUnit)); % indexes of m.pd (for every frame) before a spike initialised
        for kk = 1:length(singleUnit)
            shiftedPD = m.pd - singleUnit(kk); % subtracts spike time from pd times
            lastFrameIdx(kk) = length(shiftedPD(shiftedPD <= 0)) + opt.tsdnLatency*m.sRateHz/1000; 
                % find the indexes that are negative - then find the length
                % of just the negative - this will give you the index of
                % the frame just before the spike
        end
        
        lastFrameIdxList{ii} = lastFrameIdx;
%         stim.lastFrameIdx = lastFrameIdx;
    end
end
end

