function [] = plotHeatMap(m, s, selectUnits, tsdnLatency, sqSize, drawingMode)
%% Plot Receptive Fields
%%%% INPUT: TarTraj (above), Gaze_rest, head_offset, spikes_180 (spike times in 180Hz projector indices)
%%%% Requires: plotellipse.m

%% set up some variables

if ~strcmpi(drawingMode, 'cone') ||~strcmpi(drawingMode, 'blob') || ~strcmpi(drawingMode, 'arrows')
    drawingMode = 'arrows';
end

nLoops = m.StimGL_nloops;
Xedges = floor(min(m.angleStimXYPos(:,1))/sqSize)*sqSize:sqSize:ceil(max(m.angleStimXYPos(:,1))/sqSize)*sqSize;
Yedges = floor(min(m.angleStimXYPos(:,2))/sqSize)*sqSize:sqSize:ceil(max(m.angleStimXYPos(:,2))/sqSize)*sqSize;

for ii = selectUnits
    % heat map
    singleUnit = m.units{ii};
    singleUnit = singleUnit(m.pd(m.repeatIndex(1)) < singleUnit & singleUnit < m.pd(end));
    lastFrameIdx = ones(length(singleUnit) + 1,1);
    for jj = length(singleUnit)
        shiftedPD = m.pd(lastFrameIdx(jj):end) - singleUnit(jj);
        [~, lastFrameIdx(jj+1)] = max(shiftedPD(shiftedPD < 0));
    end
    lastFrameIdx(1) = [];
    spikeXYPos = matchedAngleStimXYPos(lastFrameIdx,:);
    spikeVel = matchedAngleStimVel(lastFrameIdx,:);
    
    figure
    hist = histogram2(spikeXYPos(:,1),spikeXYPos(:,2),Xedges,Yedges, ...
        'DisplayStyle','tile','ShowEmptyBins','on');
    xlabel('azimuth (°)')
    ylabel('elevation (°)')
    saveas(gcf, ['heatmap_unit_' num2str(ii)], 'epsc');
    
    %arrows
    figure
    lowerEdges = floor(spikeXYPos/sqSize)*sqSize;
    uniqueLowerEdges = unique(lowerEdges);
    for jj = length(uniqueLowerEdges)
        sqIndex = ismember(lowerEdges,uniqueLowerEdges(jj,:),'ropws');
        if strcmpi(drawingMode, 'arrows')
            smallQui = quiver(zeroes(size(sqIndex)), zeroes(size(sqIndex)), ...
                spikeVel(sqIndex,1)./norm(spikeVel(sqIndex,1)), ...
                spikeVel(sqIndex,2)./norm(spikeVel(sqIndex,2)), 'color',[0.2 0.2 0.2], ...
                'ShowArrowHead', 'off');
        elseif strcmpi(drawingMode, 'blob')
            spikeVecAngles = atan2d(spikeVel(sqIndex,2),spikeVel(sqIndex,1));
            leftRegion = find(spikeVecAngles < -177.5 | spikeVecAngles > 177.5);
            angleEdges = -177.5:5:177.5;
            angleMiddles = 5:5:85;
            angleHist = histogram(spikeVecAngles,angleEdges);
            normHistVals = [angleHist.Values length(leftRegion)];
            normHistVals = normHistVals/max(normHistVals);
            -normHistVals(1:17).*cos(angleMiddles);
            -normHistVals(1:17).*sin(angleMiddles);
            -normHistVals(18)
            -normHistVals(19:35).*cos(angleMiddles);
            normHistVals(19:35).*sin(angleMiddles);
            normHistVals(36)
            normHistVals(37:53).*cos(angleMiddles);
            normHistVals(37:53).*sin(angleMiddles);
            normHistVals(54)
            normHistVals(55:71).*cos(angleMiddles);
            -normHistVals(55:71).*sin(angleMiddles);
            -normHistVals(72)
            
        end
        meanQui = quiver(0, 0, 1.5*mean(spikeVel(sqIndex,1))/norm(mean(spikeVel(sqIndex,1))), ...
            1.5*mean(spikeVel(sqIndex,2))/norm(mean(spikeVel(sqIndex,2))), 'color',[1 1 1]);
        saveas(gcf, ['heatmap_unit_', num2str(ii), '_loweredges_', ...
            num2str(uniqueLowerEdges(jj,1)), '_', num2str(uniqueLowerEdges(jj,2))], 'epsc');
    end
end
%spike_tail = 5; % how many samples to plot leading up to the spike location
end