function [] = plotPolar(m, s, selectUnits, discardCorner)

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
            priorVecAngles = atan2(m.oobMatchAngVel(:,2), m.oobMatchAngVel(:,1));
            spikeVel = m.oobMatchAngVel(lastFrameIdx,:);
        else
            priorVecAngles = atan2(m.matchAngVel(:,2), m.matchAngVel(:,1));
            spikeVel = m.matchAngVel(lastFrameIdx,:);
        end
        
        spikeVecAngles = atan2(spikeVel(:,2),spikeVel(:,1));
        f = figure('visible', 'off');
        h1 = polarhistogram(priorVecAngles, 24);
        dir1=h1.Values;
        h2 = polarhistogram(spikeVecAngles, 24);
        dir2=h2.Values;
        close(f)
        
        % normalise response directions by the frequency of direction
        % presentation. Find the mean normalised binned response vector
        dir1=dir1/sum(dir1); %prior
        dir2=dir2/sum(dir2); %response
        dir3=dir2./dir1; %normalization
        dir3=dir3*100/sum(dir3);
        meanVecangle = mean(dir3,'omitnan');

        polarplot(linspace(0, 2*pi, 25),[dir3 dir3(1)]);
        hold on
        set(gcf,'color','w');
        polarplot([meanVecangle meanVecangle], [0, max(dir3)]);
        title(sprintf('unit\\_%02i',ii))
        
        saveas(gcf, ['polar_unit_' num2str(ii) '_post'], 'epsc');
        saveas(gcf, ['polar_unit_' num2str(ii) '_post'], 'fig');
    else
        warning(['Unit ' num2str(ii) ' has no spikes. A heatmap will not be plotted...']);
    end
end
end