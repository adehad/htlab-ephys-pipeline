function [] = plotHeatMap(m, s, selectUnits, tsdnLatency, sqSize, drawingMode)
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
            lastFrameIdx(kk+1) = length(shiftedPD);
        end
        lastFrameIdx(1) = [];
        spikeXYPos = m.matchedAngleStimXYPos(lastFrameIdx,:);
        spikeVel = m.matchedAngleStimVel(lastFrameIdx,:);
        figure
        histogram2(spikeXYPos(:,1),spikeXYPos(:,2),Xedges,Yedges, ...
            'DisplayStyle','tile','ShowEmptyBins','on');
        xlabel('azimuth (°)')
        ylabel('elevation (°)')
        saveas(gcf, ['heatmap_unit_' num2str(ii)], 'epsc');
        
        %arrows
        lowerEdges = floor(spikeXYPos/sqSize)*sqSize;
        uniqueLowerEdges = unique(lowerEdges,'rows');
        for jj = 1:length(uniqueLowerEdges)
            sqIndex = ismember(lowerEdges,uniqueLowerEdges(jj,:),'rows');
            sqMiddle = uniqueLowerEdges(jj,:) + sqSize/2;
            if strcmpi(drawingMode, 'arrows')
                hold on
                quiver(ones(size(find(sqIndex)))*sqMiddle(1), ones(size(find(sqIndex)))*sqMiddle(2), ...
                    2*spikeVel(sqIndex,1)./norm(spikeVel(sqIndex,1)), ...
                    2*spikeVel(sqIndex,2)./norm(spikeVel(sqIndex,2)), 'color',[0.2 0.2 0.2], ...
                    'ShowArrowHead', 'off');
            elseif strcmpi(drawingMode, 'blob')
                blobVertex = [];
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
            quiver(sqMiddle(1), sqMiddle(2), 3*mean(spikeVel(sqIndex,1))/norm(mean(spikeVel(sqIndex,1))), ...
                3*mean(spikeVel(sqIndex,2))/norm(mean(spikeVel(sqIndex,2))), 'color',[1 1 1],'ShowArrowHead', 'on');
            %saveas(gcf, ['heatmap_unit_', num2str(ii), '_loweredges_', ...
                %num2str(uniqueLowerEdges(jj,1)), '_', num2str(uniqueLowerEdges(jj,2))], 'epsc');
            hold off
        end
    else
        warning(['Unit ' num2str(ii) ' has no spikes. A heatmap will not be plotted...']);
    end
end
%spike_tail = 5; % how many samples to plot leading up to the spike location
end