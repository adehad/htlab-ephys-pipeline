function histVal = plotHeatmap(m, s, stim, selectUnits, opt, saveFig)
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
    selectUnits = 1:length(s.clusters);
end

% get bottom and left edges of each heat map bin
if strcmpi(opt.units, 'pixels')
    pos = stim.stimPos;
else
    pos = stim.angleStimPos;
end

xE = opt.sqSize*(floor(min(pos(:,1))/opt.sqSize):ceil(max(pos(:,1))/opt.sqSize));
yE = opt.sqSize*(floor(min(pos(:,2))/opt.sqSize):ceil(max(pos(:,2))/opt.sqSize));
if opt.halfOffset
    %xE = [xE - opt.sqSize/2, xE(end) + opt.sqSize/2];
    yE = [yE - opt.sqSize/2, yE(end) + opt.sqSize/2];
end
xC = xE(1:end-1)+opt.sqSize/2;
yC = yE(1:end-1)+opt.sqSize/2;

sizeMod = opt.sqSize/5; % controls sizes of arrows

%% heat map
for ii = selectUnits
    singleUnit = double(s.(sprintf('unit_%s',s.clusters(ii)))); % get sing unit spikes
    clear sqC sqMean
    if ~isempty(singleUnit)
        % get spikes for full experiment
        singleUnit = singleUnit(m.pd(1) <= singleUnit & singleUnit <= m.pd(stim.loopEndIdx(end)));
        
        % shift m.pd by the times of each spike to find out the index of last pd
        % event before the spike occurs. Introduce dragonfly neural latency
        % if wanted
        lastFrameIdx = ones(size(singleUnit));
        for kk = 1:length(singleUnit)
            shiftedPD = m.pd - singleUnit(kk);
            lastFrameIdx(kk) = length(shiftedPD(shiftedPD <= 0));% + opt.tsdnLatency*m.sRateHz/1000;
        end
        
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
        [histVal] = histcounts2(spikePos(:,1),spikePos(:,2),xE, yE);
        modHistVal = histVal/length(stim.stimLength);

        imagesc(xC,yC,modHistVal');
        axis xy % ensure y axis points up
        colorbar
        set(gcf,'color','w');
        set(gcf,'renderer','Painters');
        title(sprintf('unit\\_%s',s.clusters(ii)))
        xlabel('azimuth (°)'); ylabel('elevation (°)');
        colormap hot; c = colorbar; c.Label.String = 'Spike count per trial';
        axis equal tight
        
        % find the nearest bottom and left bin edge from each spike
        lowerE(:,1) = interp1(xE,xE,spikePos(:,1),'previous');
        lowerE(:,2) = interp1(yE,yE,spikePos(:,2),'previous');
        uniqueLE = unique(lowerE,'rows');
        
        % group spikes within same bins and find their orientation angle
        for jj = 1:size(uniqueLE,1)
            sqInd = ismember(lowerE,uniqueLE(jj,:),'rows');
            sqC(jj,:) = uniqueLE(jj,:) + opt.sqSize/2;
            sqNormVel = spikeVel(sqInd,:)./vecnorm(spikeVel(sqInd,:),2,2);
            
            if strcmpi(opt.drawingMode, 'arrows')
                hold on
                % draw small arrows from the centre of each bin for each
                % trajectory angle
                quiver(ones(size(sqNormVel))*sqC(jj,1), ones(size(sqNormVel))*sqC(jj,2), ...
                    sizeMod*sqNormVel(:,1), sizeMod*sqNormVel(:,2), ...
                    'color', 'g', 'ShowArrowHead', 'off', 'AutoScale', 'off');
            elseif strcmpi(opt.drawingMode, 'blob')
                hold on
                % draw a blob in the centre of each bin that has variable
                % radius proportional to prevalence of trajectory angle
                blobVertex = [];
                spikeVecAngles = atan2d(spikeVel(sqInd,2),spikeVel(sqInd,1));
                leftRegion = find(spikeVecAngles < -172.5 | spikeVecAngles > 172.5);
                angleE = -172.5:15:172.5;
                angleM = -165:15:180;
                
                angleHV = histcounts(spikeVecAngles,angleE);
                normHV = [angleHV length(leftRegion)];
                normHV = sizeMod*((normHV)/max(normHV) + 0.2);
                
                blobVertex(1,:) = normHV.*cosd(angleM) + sqC(jj,1);
                blobVertex(2,:) = normHV.*sind(angleM) + sqC(jj,2);
                blobVertex = [blobVertex blobVertex(:,1)];
                patch(blobVertex(1,:), blobVertex(2,:), 'g', 'EdgeColor', 'g');
            end
            
            % draw the big mean directional arrow of each bin
            hold on
            sqMean(jj,:) = mean(sqNormVel,1,'omitnan');
        end
        
        quiver(sqC(:,1), sqC(:,2), 2*sizeMod*sqMean(:,1), 2*sizeMod*sqMean(:,2), ...
                'color', 'c', 'AutoScale', 'off', 'ShowArrowHead', 'off');
            
        set(gca,'FontSize',24)
        hold off
        
        % save figures
        if saveFig == 2
            export_fig(sprintf('%s_heatmap_unit_%s_post.eps',opt.preName,s.clusters(ii)))
            saveas(gcf, join([opt.preName '_heatmap_unit_' s.clusters(ii) '_post'],''), 'fig');
        elseif saveFig
            saveas(gcf, join([opt.preName '_heatmap_unit_' s.clusters(ii) '_post'],''), 'epsc');
            saveas(gcf, join([opt.preName '_heatmap_unit_' s.clusters(ii) '_post'],''), 'fig');
        end
    else
        warning(['Unit ' s.clusters(ii) ' has no spikes. A heatmap will not be plotted...']);
    end
end
end