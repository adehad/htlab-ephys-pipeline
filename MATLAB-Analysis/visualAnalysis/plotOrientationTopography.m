function [] = plotOrientationTopography(m, s, stim, selectUnits, opt, saveFig)
%% set up some variables
if strcmpi(selectUnits, 'all')
    selectUnits = s.clusters;
end

%%
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
            spikePos = stim.oobMatchAngPos(lastFrameIdx,:);
            spikeVel = stim.oobMatchAngVel(lastFrameIdx,:);
        else
            spikePos = stim.matchAngPos(lastFrameIdx,:);
            spikeVel = stim.matchAngVel(lastFrameIdx,:);
        end
        
        % find the nearest bottom and left bin edge from each spike
        lowerEdges = floor(spikePos/opt.sqSize)*opt.sqSize;
        uniqueLowerEdges = unique(lowerEdges,'rows');
        uniqueLowerEdges(isnan(uniqueLowerEdges(:,1)),:) = [];

        if size(uniqueLowerEdges,1) > 1
        
        % group spikes within same bins and find their orientation angle
        meanVectors = zeros(size(uniqueLowerEdges));
        for jj = 1:size(uniqueLowerEdges,1)
            sqIndex = ismember(lowerEdges,uniqueLowerEdges(jj,:),'rows');
            meanVectors(jj,:) = mean(spikeVel(sqIndex,:),1);
        end
        
        % find the mean orientation angle relative to which other angles are plotted
        refAngle = atan2d(mean(meanVectors(:,2)), mean(meanVectors(:,1)));
        
        % find the relative orientation angles of each bin
        meanAngles = atan2d(meanVectors(:,2), meanVectors(:,1));
        meanAngles = meanAngles - refAngle;
        
        % lay out the mesh for the surface plot based on the centres of the bins
        [X, Y] = meshgrid(unique(uniqueLowerEdges(:,1)), unique(sort(uniqueLowerEdges(:,2))));
        X = X + opt.sqSize/2; Y = Y + opt.sqSize/2;
        Z = nan(size(X));
        
        % find the locations of each relative orientation vector on the surface
        zIdx = [];
        [~,zIdx(:,1)] = ismember(uniqueLowerEdges(:,1), unique(uniqueLowerEdges(:,1)));
        [~,zIdx(:,2)] = ismember(uniqueLowerEdges(:,2), unique(sort(uniqueLowerEdges(:,2))));
        for jj = 1:size(zIdx,1)
            Z(zIdx(jj,2),zIdx(jj,1)) = meanAngles(jj);
        end
        
        % plot the surface
        figure
        set(gcf,'color','w');
        surf(X,Y,Z);
        colorbar
        title(sprintf('unit\\_%02i',ii))
        
        % save figures
        if saveFig
            saveas(gcf, [opt.preName '_orientationTopo_unit_' num2str(ii) '_post'], 'epsc');
            saveas(gcf, [opt.preName '_orientationTopo_unit_' num2str(ii) '_post'], 'fig');
        end
        else
            warning(['Unit ' num2str(ii) ' has only one filled bin. Topography will not be plotted...']);
        end
    else
        warning(['Unit ' num2str(ii) ' has no spikes. A heatmap will not be plotted...']);
    end
end
end