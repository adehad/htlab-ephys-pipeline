function [] = plotHeatMap(m, s, selectUnits, tsdnLatency, sqSize, drawingMode, discardCorner)
%% set up some variables

if ~strcmpi(drawingMode, 'cone') && ~strcmpi(drawingMode, 'blob') && ~strcmpi(drawingMode, 'arrows')
    drawingMode = 'arrows';
end

if strcmpi(selectUnits, 'all')
    selectUnits = s.clusters;
end

Xedges = floor(min(m.angleStimXYPos(:,1))/sqSize)*sqSize:sqSize:ceil(max(m.angleStimXYPos(:,1))/sqSize)*sqSize;
Yedges = floor(min(m.angleStimXYPos(:,2))/sqSize)*sqSize:sqSize:ceil(max(m.angleStimXYPos(:,2))/sqSize)*sqSize;

for ii = selectUnits
    % heat map
    singleUnit = double(s.(sprintf('unit_%02i',ii)));
    if ~isempty(singleUnit)
        singleUnit = singleUnit(m.pd(1) <= singleUnit & singleUnit <= m.pd(m.repeatIndex(end)));
        lastFrameIdx = ones(length(singleUnit) + 1,1);
        for kk = 1:length(singleUnit)
            shiftedPD = m.pd(lastFrameIdx(kk):end) - singleUnit(kk);
            shiftedPD = shiftedPD(shiftedPD <= 0);
            lastFrameIdx(kk+1) = length(shiftedPD) + tsdnLatency;
        end
        lastFrameIdx(1) = [];
        
        mASPos = m.matchedAngleStimXYPos;
        mASV = m.matchedAngleStimVel;
        if discardCorner
            ooB = m.outOfBoundsIdx+1;
            ooB(find(ooB == 2)) = nan;
            mASPos = mASPos.*ooB;
            mASV = mASV.*ooB;
        end
        spikeXYPos = mASPos(lastFrameIdx,:);
        spikeVel = mASV(lastFrameIdx,:);
        
        figure
        histogram2(spikeXYPos(:,1),spikeXYPos(:,2),Xedges,Yedges, ...
            'DisplayStyle','tile','ShowEmptyBins','on');
        xlabel('azimuth (°)')
        ylabel('elevation (°)')
        colormap hot
        colorbar
        axis equal
        %saveas(gcf, ['heatmap_unit_' num2str(ii) '_pre'], 'epsc');

        lowerEdges = floor(spikeXYPos/sqSize)*sqSize;
        uniqueLowerEdges = unique(lowerEdges,'rows');
        for jj = 1:size(uniqueLowerEdges,1)
            sqIndex = ismember(lowerEdges,uniqueLowerEdges(jj,:),'rows');
            sqMiddle = uniqueLowerEdges(jj,:) + sqSize/2;
            if strcmpi(drawingMode, 'arrows')
                hold on
                quiver(ones(size(find(sqIndex)))*sqMiddle(1), ones(size(find(sqIndex)))*sqMiddle(2), ...
                    2*spikeVel(sqIndex,1)./vecnorm(spikeVel(sqIndex,:),2,2), ...
                    2*spikeVel(sqIndex,2)./vecnorm(spikeVel(sqIndex,:),2,2), ...
                    'color', 'g', 'ShowArrowHead', 'off', 'AutoScale', 'off');
            elseif strcmpi(drawingMode, 'blob')
                hold on
                blobVertex = [];
                spikeVecAngles = atan2d(spikeVel(sqIndex,2),spikeVel(sqIndex,1));
                leftRegion = find(spikeVecAngles < -172.5 | spikeVecAngles > 172.5);
                angleEdges = -172.5:15:172.5;
                angleMiddles = -165:15:180;
                
                angleHist = histcounts(spikeVecAngles,angleEdges);
                normHistVals = [angleHist length(leftRegion)];
                normHistVals = 2*((normHistVals)/max(normHistVals) + 0.2);
                
                blobVertex(1,:) = normHistVals.*cosd(angleMiddles) + sqMiddle(1);
                blobVertex(2,:) = normHistVals.*sind(angleMiddles) + sqMiddle(2);
                blobVertex = [blobVertex blobVertex(:,1)];
                patch(blobVertex(1,:), blobVertex(2,:), 'g', 'EdgeColor', 'g');
            end
            meanVec = mean(spikeVel(sqIndex,:),1);
            quiver(sqMiddle(1), sqMiddle(2), 4*meanVec(1)/norm(meanVec), ...
                4*meanVec(2)/norm(meanVec), ...
                'color', 'c', 'MaxHeadSize', 7, 'AutoScale', 'off');
            %saveas(gcf, ['heatmap_unit_', num2str(ii), '_loweredges_', ...
                %num2str(uniqueLowerEdges(jj,1)), '_', num2str(uniqueLowerEdges(jj,2))], 'epsc');
            hold off
        end
        saveas(gcf, ['heatmap_unit_' num2str(ii) '_post'], 'epsc');
        saveas(gcf, ['heatmap_unit_' num2str(ii) '_post'], 'fig');
    else
        warning(['Unit ' num2str(ii) ' has no spikes. A heatmap will not be plotted...']);
    end
end
end