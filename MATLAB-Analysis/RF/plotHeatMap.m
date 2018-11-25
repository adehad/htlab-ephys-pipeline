function [] = plotHeatMap(m, s, selectUnits, tsdnLatency, sqSize, drawingMode)
%% set up some variables

if ~strcmpi(drawingMode, 'cone') ||~strcmpi(drawingMode, 'blob') || ~strcmpi(drawingMode, 'arrows')
    drawingMode = 'arrows';
end

if strcmpi(selectUnits, 'all')
    selectUnits = s.cluster_groups.cluster_id;
end

Xedges = floor(min(m.angleStimXYPos(:,1))/sqSize)*sqSize:sqSize:ceil(max(m.angleStimXYPos(:,1))/sqSize)*sqSize;
Yedges = floor(min(m.angleStimXYPos(:,2))/sqSize)*sqSize:sqSize:ceil(max(m.angleStimXYPos(:,2))/sqSize)*sqSize;

for ii = selectUnits
    % heat map
    singleUnit = double(s.(sprintf('unit_%02i',ii)));
    singleUnit = singleUnit(m.pd(1) <= singleUnit & singleUnit <= m.pd(m.repeatIndex(end)));
    lastFrameIdx = ones(length(singleUnit) + 1,1);
    for kk = 1:length(singleUnit)
        shiftedPD = m.pd(lastFrameIdx(kk):end) - singleUnit(kk);
        shiftedPD = shiftedPD(shiftedPD <= 0);
        lastFrameIdx(kk+1) = length(shiftedPD);
    end
    lastFrameIdx(1) = [];
    spikeXYPos = matchedAngleStimXYPos(lastFrameIdx,:);
    spikeVel = matchedAngleStimVel(lastFrameIdx,:);
    figure
    histogram2(spikeXYPos(:,1),spikeXYPos(:,2),Xedges,Yedges, ...
        'DisplayStyle','tile','ShowEmptyBins','on');
    xlabel('azimuth (°)')
    ylabel('elevation (°)')
    saveas(gcf, ['heatmap_unit_' num2str(ii)], 'epsc');
    
    %arrows
    lowerEdges = floor(spikeXYPos/sqSize)*sqSize;
    uniqueLowerEdges = unique(lowerEdges);
    for jj = length(uniqueLowerEdges)
        figure
        sqIndex = ismember(lowerEdges,uniqueLowerEdges(jj,:),'rows');
        if strcmpi(drawingMode, 'arrows')
            hold on
            quiver(zeroes(size(sqIndex)), zeroes(size(sqIndex)), ...
                spikeVel(sqIndex,1)./norm(spikeVel(sqIndex,1)), ...
                spikeVel(sqIndex,2)./norm(spikeVel(sqIndex,2)), 'color',[0.2 0.2 0.2], ...
                'ShowArrowHead', 'off');
        elseif strcmpi(drawingMode, 'blob')
            spikeVecAngles = atan2d(spikeVel(sqIndex,2),spikeVel(sqIndex,1));
            leftRegion = find(spikeVecAngles < -177.5 | spikeVecAngles > 177.5);
            angleEdges = -177.5:5:177.5;
            angleMiddles = -175:5:180;
            
            angleHist = histogram(spikeVecAngles,angleEdges);
            normHistVals = [angleHist.Values length(leftRegion)];
            normHistVals = normHistVals/max(normHistVals);
            
            blobVertex(1,:) = normHistVals.*cos(angleMiddles);
            blobVertex(2,:) = normHistVals.*sin(angleMiddles);
            hold on
            for kk = length(blobVertex)-1
                plot(blobVertex(1,[kk kk+1]),blobVertex(2,[kk kk+1]), 'color',[0.5 0.5 0.5]);
            end
            plot(blobVertex(1,[end 1]),blobVertex(2,[end 1]), 'color',[0.5 0.5 0.5]);
        end
        quiver(0, 0, 1.5*mean(spikeVel(sqIndex,1))/norm(mean(spikeVel(sqIndex,1))), ...
            1.5*mean(spikeVel(sqIndex,2))/norm(mean(spikeVel(sqIndex,2))), 'color',[1 1 1]);
        saveas(gcf, ['heatmap_unit_', num2str(ii), '_loweredges_', ...
            num2str(uniqueLowerEdges(jj,1)), '_', num2str(uniqueLowerEdges(jj,2))], 'epsc');
        hold off
    end
end
%spike_tail = 5; % how many samples to plot leading up to the spike location
end