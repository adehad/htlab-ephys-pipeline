function [events,waves,xi] = get_events(x,ch,f,win,msec,thr,bgoff,recenter)
%
% given a time series, find events > than a threshold
%
%IN:  x  = time series; could be matrix
%    ch  = channel of interest
%     f  = denoising filter
%   win  = +- window to extract around events (in sample pts!)
%   msec = pts/msec
%    thr = event detection threshold for *filtered* bgoff subtracted waveforms
%  bgoff = background subtract before threshold (1=true)
%  recenter = recenter detected events to event peak (1=true)
%
%OUT: events = vector list of event times
%     waves = matrix of detected event waveforms of +-win + 1 in length
%         xi = filtered time series for x
%
% AL, janelia 9/2010
%
% [events,waves,xi] = get_events(x,ch,f,win,msec,thr,bgoff,recenter)
% 
%Requires: splitconv.m & tconv.m

% misc params
flen = length(f);
tmax = size(x,2);

% get events
xi = splitconv(x(ch,:),f);
if bgoff
    xi = xi - mean(xi(flen:end-flen));
end;

if thr > 0
% events = find(diff(xi > thr) == 1) + 1;
events = find(diff(xi > thr) > 0) + 1; %a spike event is when the signal crosses the thr value
else
events = find(diff(xi > thr) < 0) + 1; %a spike event is when the signal crosses the thr value    
end

% get waves
waves = [];
preAlloc = 1e4; % 10k pre allocated rows, makes it a bit speedy
waves = zeros(preAlloc,2*win+1); 
recenter = 3; %set this to 3 here for now
if recenter ~= 0 %if recentering is required, we run the wave cropping loop twice
    tloops = 2
%     disp('recentering!!')    
else tloops = 1 ;
end;

for j = 1:length(events)
    event_j = events(j);
    for k = 1:tloops
        if event_j-win < 1
            prew = 1;
            preb = win - event_j + 1;
        else prew = event_j - win;
            preb = 0;
        end;
        if event_j+win > tmax
            postw = tmax;
            postb = win - (tmax - event_j);
        else postw = event_j+win;
            postb = 0;
        end;
        wavetmp = [zeros(1,preb) x(ch,prew:postw) zeros(1,postb)];
        if bgoff
            wavetmp = wavetmp - mean(wavetmp(win-round(0.75*msec):win-round(0.5*msec)));
        end;
        if k==1;
            switch recenter
                case 1 %center at absolute maximum
                    tnew = find(abs(wavetmp(round(win-0.5*msec):round(win+0.5*msec))) == max(abs(wavetmp(round(win-0.5*msec):round(win+0.5*msec))))) + win-0.5*msec - 1;
                case 2 %center at the midpoint of the triggered positive peak
                    tnew = find(wavetmp(round(win-0.5*msec):round(win+0.5*msec)) > thr,1) + win-0.5*msec - 1;
                case 3 %center at the negative peak
                    tnew = find(wavetmp(win:round(win+0.25*msec)) == min(wavetmp(win:round(win+0.25*msec)))) + win;
%                     disp('neg peak')
                otherwise  tnew = [];
            end;
        end
        if ~isempty(tnew) %correction is required
            event_j = tnew(1) - (win+1) + event_j - 1;
%         else
%             clf;
%             plot(wavetmp)
%             hold on
%             pause
%             sprintf('tnew null!!')
        end;
    end;
    events(j) = event_j;
%     waves = [waves; wavetmp];
    waves(j,:) = wavetmp;
end;

waves = waves(1:j,:);
