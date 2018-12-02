function [] = plotPolarPlot(m, s, selectUnits, discardCorner)

if strcmpi(selectUnits, 'all')
    selectUnits = s.clusters;
end

for ii = selectUnits
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
        
        if discardCorner
            spikeVel = m.oobMatchAngVel(lastFrameIdx,:);
        else
            spikeVel = m.matchAngVel(lastFrameIdx,:);
        end
        
        figure
        set(gcf,'color','w');
        
        spikeVecAngles = atan2d(spikeVel(:,2),spikeVel(:,1));
        leftRegion = find(spikeVecAngles < -172.5 | spikeVecAngles > 172.5);
        angleEdges = -172.5:15:172.5;
        angleMiddles = deg2rad([195:15:360, 15:15:195]);
        
        angleHist = histcounts(spikeVecAngles,angleEdges);
        normHistVals = [angleHist length(leftRegion)];
        normHistVals = 4*((normHistVals)/max(normHistVals) + 0.2);
        normHistVals = [normHistVals normHistVals(1)];
        polarplot(angleMiddles, normHistVals);
        
        %TODO: this mean is done on -180 180
        hold on
        meanVecangle = mean(spikeVecAngles,'omitnan');
        if meanVecangle < 0
            meanVecangle = meanVecangle + 360;
        end
        polarplot(deg2rad([meanVecangle meanVecangle]), [0, max(normHistVals)]);
        
        title(sprintf('unit\\_%02i',ii))
        saveas(gcf, ['polar_unit_' num2str(ii) '_post'], 'epsc');
        saveas(gcf, ['polar_unit_' num2str(ii) '_post'], 'fig');
    else
        warning(['Unit ' num2str(ii) ' has no spikes. A heatmap will not be plotted...']);
    end
end
end