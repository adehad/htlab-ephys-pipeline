function [] = plotPolar(m, s, stim, selectUnits, opt, saveFig)

if strcmpi(selectUnits, 'all')
    selectUnits = 1:length(s.clusters);
end

for ii = selectUnits
    singleUnit = double(s.(sprintf('unit_%s',s.clusters(ii))));
    if ~isempty(singleUnit)
        singleUnit = singleUnit(m.pd(1) <= singleUnit & singleUnit <= m.pd(stim.repeatIndex(end)));
        lastFrameIdx = ones(length(singleUnit) + 1,1);
        for kk = 1:length(singleUnit)
            shiftedPD = m.pd(lastFrameIdx(kk):end) - singleUnit(kk);
            shiftedPD = shiftedPD(shiftedPD <= 0);
            lastFrameIdx(kk+1) = length(shiftedPD);
        end
        lastFrameIdx(1) = [];
        
        if opt.discardCorner
            priorVecAngles = atan2(stim.oobMatchAngVel(:,2), stim.oobMatchAngVel(:,1));
            spikeVel = stim.oobMatchAngVel(lastFrameIdx,:);
        else
            priorVecAngles = atan2(stim.matchAngVel(:,2), stim.matchAngVel(:,1));
            spikeVel = stim.matchAngVel(lastFrameIdx,:);
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
        title(sprintf('unit\\_%s',s.clusters(ii)))
        hold off
        
        if saveFig == 2
            export_fig(sprintf('%s_polar_unit_%s.eps',opt.preName,num2str(ii)))
            saveas(gcf, [opt.preName '_polar_unit_' s.clusters(ii)], 'fig');
        elseif saveFig
            saveas(gcf, [opt.preName '_polar_unit_' s.clusters(ii)], 'epsc');
            saveas(gcf, [opt.preName '_polar_unit_' s.clusters(ii)], 'fig');
        end
    else
        warning(['Unit ' s.clusters(ii) ' has no spikes. A heatmap will not be plotted...']);
    end
end
end