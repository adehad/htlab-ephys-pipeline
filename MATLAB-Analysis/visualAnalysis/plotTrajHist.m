function [n, q, b] = plotTrajHist(m, s, stim, selectUnits, opt, saveFig)
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

%% heat map
u = 0;
for ii = selectUnits
    singleUnit = double(s.(sprintf('unit_%s',s.clusters(ii)))); % get sing unit spikes
    if ~isempty(singleUnit)
        % get spikes for full experiment
        singleUnit = singleUnit(m.pd(1) <= singleUnit & singleUnit <= m.pd(stim.repeatIndex(end)));
        
        % shift m.pd by the times of each spike to find out the index of last pd
        % event before the spike occurs. Introduce dragonfly neural latency
        % if wanted
        lastFrameIdx = ones(size(singleUnit));
        for kk = 1:length(singleUnit)
            shiftedPD = m.pd - singleUnit(kk);
            shiftedPD = shiftedPD(shiftedPD <= 0);
            lastFrameIdx(kk) = length(shiftedPD(shiftedPD <= 0)) + opt.tsdnLatency*m.sRateHz/1000;
        end
        
        % discard data in the corner of the screen where the target goes to
        % when the trajectory goes out of bounds, if wanted
        spikeTraj = stim.pdTraj(lastFrameIdx,:);
        
        e = linspace(0,1,50);
        c = e(1:end-1)+e(2)/2;
        for kk = 0:2:6
            figure
            hold on
            for jj = [1 9 17 25 33 41]+kk
                u = u + 1;
                t = find(spikeTraj(:,1) == jj);
                q{u} = spikeTraj(t,2);
                if jj == 9+kk || jj == 25+kk || jj == 41+kk
                    q{u} = 1-q{u};
                end
                [n(u,:), ~, b{u}] = histcounts(q{u},e);
                plot(c,n(u,:)/length(stim.stimLength))
            end
            set(gcf,'color','w');
            set(gcf,'renderer','Painters');
            title(sprintf('unit\\_%s\\_dir\\_%s',s.clusters(ii),num2str(kk)))
            xlabel('Normalised distance along trajectory line'); ylabel('Spike count per trial');
            legend()
            hold off
        end
    end
end
end