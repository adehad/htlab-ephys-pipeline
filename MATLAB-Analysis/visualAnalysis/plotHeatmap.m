function [] = plotHeatmap(m, s, stim, selectUnits, opt, saveFig)
%% set up some variables

%TODO: SCALE MEAN ARROWS

% choose spike direction drawing type
if ~strcmpi(opt.drawingMode, 'blob') && ~strcmpi(opt.drawingMode, 'arrows')
    opt.drawingMode = 'none';
end
% choose whether heat map is on screen pixels or angular
if ~strcmpi(opt.units, 'pixels') && ~strcmpi(opt.units, 'angle')
    opt.units = 'angle';
end
if strcmpi(selectUnits, 'all')
    selectUnits = s.clusters;
end

% get bottom and left edges of each heat map bin
if strcmpi(opt.units, 'pixels')
    Xedges = floor(min(stim.stimPos(:,1)-1)/opt.sqSize)*opt.sqSize:opt.sqSize:...
        ceil(max(stim.stimPos(:,1)+1)/opt.sqSize)*opt.sqSize;
    Yedges = floor(min(stim.stimPos(:,2)-1)/opt.sqSize)*opt.sqSize:opt.sqSize:...
        ceil(max(stim.stimPos(:,2)+1)/opt.sqSize)*opt.sqSize;
else
    Xedges = floor(min(stim.angleStimPos(:,1))/opt.sqSize)*opt.sqSize:opt.sqSize:...
        ceil(max(stim.angleStimPos(:,1))/opt.sqSize)*opt.sqSize;
    Yedges = floor(min(stim.angleStimPos(:,2))/opt.sqSize)*opt.sqSize:opt.sqSize:...
        ceil(max(stim.angleStimPos(:,2))/opt.sqSize)*opt.sqSize;
end

sizeMod = opt.sqSize/5; % controls sizes of arrows

%% heat map
for ii = selectUnits
    singleUnit = double(s.(sprintf('unit_%02i',ii))); % get sing unit spikes
    if ~isempty(singleUnit)
        % get spikes for full experiment
        singleUnit = singleUnit(m.pd(1) <= singleUnit & singleUnit <= m.pd(stim.repeatIndex(end)));
        
        % shift m.pd by the times of each spike to find out the last pd
        % event before the spike occurs. Introduce dragonfly neural latency
        % if wanted
        lastFrameIdx = ones(length(singleUnit) + 1,1);
        for kk = 1:length(singleUnit)
            shiftedPD = m.pd(lastFrameIdx(kk):end) - singleUnit(kk);
            shiftedPD = shiftedPD(shiftedPD <= 0);
            lastFrameIdx(kk+1) = length(shiftedPD) + opt.tsdnLatency;
        end
        lastFrameIdx(1) = [];
        
        % discard data in the corner of the screen where the target goes to
        % when the trajectory goes out of bounds, if wanted
        if opt.discardCorner
            if strcmpi(opt.units, 'pixels')
                spikePos = stim.oobMatchPos(lastFrameIdx,:);
                spikeVel = stim.oobMatchVel(lastFrameIdx,:);
            else
                spikePos = stim.oobMatchAngPos(lastFrameIdx,:);
                spikeVel = stim.oobMatchAngVel(lastFrameIdx,:);
            end
        else
            if strcmpi(opt.units, 'pixels')
                spikePos = stim.matchPos(lastFrameIdx,:);
                spikeVel = stim.matchVel(lastFrameIdx,:);
            else
                spikePos = stim.matchAngPos(lastFrameIdx,:);
                spikeVel = stim.matchAngVel(lastFrameIdx,:);
            end
        end
        
        % plot histogram of spike positions
        figure
        histogram2(spikePos(:,1),spikePos(:,2),Xedges,Yedges, ...
            'DisplayStyle','tile','ShowEmptyBins','on');
        set(gcf,'color','w');
        title(sprintf('unit\\_%02i',ii))
        xlabel('azimuth (°)'); ylabel('elevation (°)');
        colormap hot; c = colorbar; c.Label.String = 'Spike count';
        axis equal
        
%         if saveFig
%             saveas(gcf, [opt.preName '_heatmap_unit_' num2str(ii) '_pre'], 'epsc');
%         end
        
        % find the nearest bottom and left bin edge from each spike
        lowerEdges = floor(spikePos/opt.sqSize)*opt.sqSize;
        uniqueLowerEdges = unique(lowerEdges,'rows');
        
        % group spikes within same bins and find their orientation angle
        for jj = 1:size(uniqueLowerEdges,1)
            sqIndex = ismember(lowerEdges,uniqueLowerEdges(jj,:),'rows');
            sqMiddle = uniqueLowerEdges(jj,:) + opt.sqSize/2;
            
            if strcmpi(opt.drawingMode, 'arrows')
                hold on
                % draw small arrows from the centre of each bin for each
                % trajectory angle
                quiver(ones(size(find(sqIndex)))*sqMiddle(1), ones(size(find(sqIndex)))*sqMiddle(2), ...
                    2*spikeVel(sqIndex,1)./vecnorm(spikeVel(sqIndex,:),2,2), ...
                    2*spikeVel(sqIndex,2)./vecnorm(spikeVel(sqIndex,:),2,2), ...
                    'color', 'g', 'ShowArrowHead', 'off', 'AutoScale', 'off');
            elseif strcmpi(opt.drawingMode, 'blob')
                hold on
                % draw a blob in the centre of each bin that has variable
                % radius proportional to prevalence of trajectory angle
                blobVertex = [];
                spikeVecAngles = atan2d(spikeVel(sqIndex,2),spikeVel(sqIndex,1));
                leftRegion = find(spikeVecAngles < -172.5 | spikeVecAngles > 172.5);
                angleEdges = -172.5:15:172.5;
                angleMiddles = -165:15:180;
                
                angleHist = histcounts(spikeVecAngles,angleEdges);
                normHistVals = [angleHist length(leftRegion)];
                normHistVals = sizeMod*((normHistVals)/max(normHistVals) + 0.2);
                
                blobVertex(1,:) = normHistVals.*cosd(angleMiddles) + sqMiddle(1);
                blobVertex(2,:) = normHistVals.*sind(angleMiddles) + sqMiddle(2);
                blobVertex = [blobVertex blobVertex(:,1)];
                patch(blobVertex(1,:), blobVertex(2,:), 'g', 'EdgeColor', 'g');
            end
            
            % draw the big mean directional arrow of each bin
            hold on
            meanVec = mean(spikeVel(sqIndex,:),1,'omitnan');
            quiver(sqMiddle(1), sqMiddle(2), 2*sizeMod*meanVec(1)/norm(meanVec), ...
                2*sizeMod*meanVec(2)/norm(meanVec), ...
                'color', 'c', 'MaxHeadSize', 7, 'AutoScale', 'off');
            
%             if saveFig
%                 saveas(gcf, [opt.preName '_heatmap_unit_', num2str(ii), '_loweredges_', ...
%                     num2str(uniqueLowerEdges(jj,1)), '_', num2str(uniqueLowerEdges(jj,2))], 'epsc');
%             end
        end
        set(gca,'FontSize',24)
        hold off
        % save figures
        if saveFig == 2
            export_fig(sprintf('%s_heatmap_unit_%s_post.eps',opt.preName,num2str(ii)))
            saveas(gcf, [opt.preName '_heatmap_unit_' num2str(ii) '_post'], 'fig');
        elseif saveFig
            saveas(gcf, [opt.preName '_heatmap_unit_' num2str(ii) '_post'], 'epsc');
            saveas(gcf, [opt.preName '_heatmap_unit_' num2str(ii) '_post'], 'fig');
        end
        
    else
        warning(['Unit ' num2str(ii) ' has no spikes. A heatmap will not be plotted...']);
    end
end
end