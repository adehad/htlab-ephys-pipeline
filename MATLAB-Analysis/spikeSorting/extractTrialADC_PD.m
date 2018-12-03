function [m] = extractTrialADC_PD(rawBinary, m, outputFileName )
% 
% Example Usage
%% Check Required Functions
% Appended to the end of this fuction
%       get_ch_thr   get_events(get_ch_thr)   splitconv(get_events)
%       tconv(splitconv)
%% Default Settings
% If the variables do not exist in the workspace, or are empty, set defaults      
 
pdwin = 4;      
block = 1e6;    
esize = 16;     
%% set up some parameters in structure m
    % Path can also be added to filenames below

% Get binary filesize
fid = fopen(rawBinary,'r');
if fid > 0
    fseek(fid,0,'eof');
    m.fileSizeBytes = ftell(fid);
    fseek(fid,0,'bof');
else
    fclose(fid);
    error(['Could not open binary file: ', rawBinary])
end

m.pdwin     = pdwin*m.msec; % window to crop around pd events
esize_msec  = esize*m.msec; % in units of sample pts
esize_bytes = esize_msec*m.nChans*m.dbytes; % in units of bytes
dlen        = m.fileSizeBytes/(m.dbytes*m.nChans);
trueblock   = block - 2*esize_msec;     % Block reading code has fseek the moves back by 2*esize_msec

%% extract waves according to spike times to confirm that the sorting is correct

% denoising filters
m.pd_flen = 200;
m.pd_cut = m.fps/2;
m.pd_f = fir1(m.pd_flen,m.pd_cut/(m.sRateHz/2),'low'); % photodiode


% set event thresholds
try         % check if m.pdthr is set
    m.pdthr;
catch       % if not set, use get_ch_thr - to set the channel threshold
    m.pdthr = get_ch_thr(rawBinary,m.nChans,m.pdch,5e5,'pd',m.pd_f,3,0,0);
end
 
x = [];        
pd = [];
pdRise = [];
pdFall = [];
pwaves = [];

% slong_srate = 1000; %Hz
% slong_drate = round(m.sRateHz/slong_srate);
% m.slong_srate = m.sRateHz/slong_drate;
% slong_block = trueblock/slong_drate;
% m.slong = int16(zeros(length(m.ech),slong_block*floor(dlen/trueblock)));

blocklen_k = 1e6;
fid = fopen(rawBinary,'r');
lastblock = floor(dlen/trueblock) + 1;

% figure(111); hold off; clf
for k = 1:lastblock  % go to the last block
    
    if k == lastblock
        blocklen_k = inf;
    else
        blocklen_k = block;
    end
    x = fread(fid,[m.nChans blocklen_k],'int16'); % load one block at a time
       
    
    % get pd data
    [pdk,pdrk,pdfk,pwavetmp] = get_events(x,m.pdch,m.pd_f,m.pdwin,m.msec,m.pdthr,0,0);
    truepd = find((pdk > esize_msec) & (pdk < size(x,2) - esize_msec));
    pwaves = [pwaves; pwavetmp(truepd,:)];
    pd = [pd (pdk(truepd)+(k-1)*(trueblock))];
    
    truepdr = find((pdrk > esize_msec) & (pdrk < size(x,2) - esize_msec));
    pdRise = [pdRise (pdrk(truepdr)+(k-1)*(trueblock))];
    truepdf = find((pdfk > esize_msec) & (pdfk < size(x,2) - esize_msec));
    pdFall = [pdFall (pdfk(truepdf)+(k-1)*(trueblock))];
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % get PD data different approach (?)
%     
%         % Filter photodiode channels - done in get_events
% %     for cc=m.pdch(1):m.pdch(end) 
% %         x(cc,:) = splitconv(x(cc,:),m.pd_f);
% %     end
%     
%     
%     % Establish b_Offset - the starting element of the block - so we can search trial_spikes properly 
%     if k < lastblock
%         b_Offset = (k-1)*blocklen_k;
%         fseek_Offset = (k-1)*(+2*esize_msec); 
%         b_Offset = b_Offset - fseek_Offset;
%     else % Comment this else part out if you do not want to search for waves in the final block
%         blocklen_k = block;
%         b_Offset = (k-1)*blocklen_k;
%         fseek_Offset = (k-1)*(+2*esize_msec); 
%         b_Offset = b_Offset - fseek_Offset;
%         blocklen_k = size(x,2);
%     end
% 
%     [pdk,pwavetmp] = get_events(x,m.pdch,m.pd_f,m.pdwin,m.msec,m.pdthr,0,0); % get indexes for pd events
%     
%     % Ensure we get the most number of extractable waveforms - by making sure
%     % event index has enough space to extract waveform
%     if k ==1
%         truepd = find( (pdk > (m.pdwin+1) ) & (pdk < (size(x,2)-esize_msec)) );
%             % First block should only look after m.spikewin to prevent
%             % negative spikes, -esize_msec is because blocks overlap
%     elseif k==lastblock
%         truepd = find( (pdk > esize_msec) & (pdk < (size(x,2)-m.spikewin-1)) );       
%             % Last block should look until m.spikewin before the end to
%             % prevent positive spikes
%     else
%         truepd = find((pdk > esize_msec) & (pdk < size(x,2) - esize_msec));
%                 % trial_spikes equivalent to pdk (event times)
%                 % +/-(esize_msec)... is because blocks overlap due to fseek - we don't want to record same spike multiple times
%     end
%     pwaves = [pwaves; pwavetmp(truepd,:)];
%     pd = [pd (pdk(truepd)+(k-1)*(trueblock))];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    fprintf('Data block %02d of %d; %d pd events found \n' ...
            ,k,floor(dlen/block)+1, length(pdk))
                                
    % Shift back - creates overlapping blocks
    fseek(fid,-2*esize_bytes,'cof');
        
%     figure(111); plot(squeeze(waves(:,:,(size(waves,3)-size(tempWaves,3)+1):50:end))); title(['Data Blocks 1 to ', num2str(k)]); hold on; drawnow; 
end
fclose(fid);
disp('-----------------------------------------------------------')


%% Save PD data to mat file

m.pd = pd;
m.pdRise = pdRise;
m.pdFall = pdFall;

%% Save to .mat File

% clear up variables before saving to keep it neat
if ~isempty(outputFileName)
    fprintf('Saving function output ... ');
    clear fid ans pwavetmp x blocklen_k k cc fseek_Offset b_Offset lastblock
    save(outputFileName)
    fprintf('done! \n');
end
end


% Appended Functions
function [events, eRising, eFalling, waves,xi] = get_events(x,ch,f,win,msec,thr,bgoff,recenter)
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
events = find(abs(diff(xi > thr)) > 0) + 1; %a spike event is when the signal crosses the thr value
eRising = find(diff(xi > thr) > 0) + 1;
eFalling = events(~ismember(events, eRising));
else
events = find(diff(xi > thr) < 0) + 1; %a spike event is when the signal crosses the thr value    
end

% get waves
waves = [];
preAlloc = 1e4; % 10k pre allocated rows, makes it a bit speedy
waves = zeros(preAlloc,2*win+1); 
recenter = 3; %set this to 3 here for now
if recenter ~= 0 %if recentering is required, we run the wave cropping loop twice
    tloops = 2;
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
end

function thr = get_ch_thr(fname,nch,tch,ndata,chname,f,thr,bgoff,recenter)
%
%allow user to set event detection threshold for a channel
%
%IN: fname = file to load
%    nch = total ch in file
%    tch = channel to set threshold
%    ndata = data load block size (pts)
%    chname = string name of ch
%    f = denoising filter for ch
%    thr = initial thr for ch
%    bgoff = subtract bg (1=yes, 0=no)
%
%OUT: thr = final thr for ch
%
%Requires: get_events.m

flen = length(f);
fid = fopen(fname,'r');
cmdlist = sprintf('channel: %s\n\n(s)et new threshold\n(l)oad next data block\n(a)ccept settings\nset (b)lock size',chname);
cmd = 'l';
while cmd ~= 'a'
    switch cmd
        case 'l'
            x = fread(fid,[nch ndata],'int16');
            if bgoff
                xmu = mean(x(tch,:));
            else xmu = 0;
            end;
            [events,waves,xi] = get_events(x,tch,f,10,10,thr,bgoff,recenter);
            clf;
            hold on
            plot(x(tch,:)-xmu,'b')
            plot(xi,'r')
            plot(events,xi(events),'m*');
            plot([0 ndata],[1 1]*thr,'k:')
            zoom on;
            title(sprintf('%s\nfound %d events',cmdlist,length(events)));
        case 's'
            thr = input(sprintf('Enter new threshold (old= %d): ',thr));
            [events,waves,xi] = get_events(x,tch,f,10,10,thr,bgoff,recenter);
            clf;
            hold on
            plot(x(tch,:)-xmu,'b')
            plot(xi,'r')
            plot(events,xi(events),'m*');
            plot([0 ndata],[1 1]*thr,'k:')
            zoom on;
            title(sprintf('%s\nfound %d events',cmdlist,length(events)));
        case 'b'
            ndata = input(sprintf('Enter new data load block size (old= %d): ',ndata));
    end;
    %     [x1,y1,cmd] = ginput2(1);
    cmd = input('next command: ','s');
end;
fclose(fid);
end

function [dnew,fmax] = tconv(data,kernel,cmid)
%
% conv data with kernel and compensate for time-shift
% induced by convolution
%
% length(dnew) = length(data)
%
% assume kernel center is at max unless cmid=1
%
%IN: data = time series to be convolved
%   kernel = convolution kernel
%   cmid = assume middle of kernel is center (ie for derivative)
%
%OUT: dnew = convolved, time-shifted data
%     fmax = center of kernel
%
% AL, Caltech, 8/00
%
% [dnew,fmax] = tconv(data,kernel,cmid)
%

if ~exist('cmid','var')  || isempty(cmid),  cmid = 0; end;
dnew = conv(data,kernel);
if (cmid == 0)
    fmax = find(kernel == max(kernel),1);
else fmax = round(length(kernel)/2);
end;
dnew = dnew(fmax:fmax+length(data)-1);

% recenter even-length (non-symmetric) filter with cmid centering)
if (mod(length(kernel),2) == 0) & cmid == 1
    dnewi = interp(dnew,2);
    dnew = dnewi(2:2:end);
end;
end

function y = splitconv(d,f)
%
% convolution with minimized edge artifacts
% assumes length(d) >> length(f)
%
%IN: d = time series (vector)
%    f = kernel
%
%OUT: y = filtered d
%
% AL, janelia; 11/09
%
%y = splitconv(d,f)
%


dmid = round(size(d,1)/2);
da = d - d(1);
db = d - d(end);
daf = tconv(da,f) + d(1);
dbf = tconv(db,f) + d(end);
daf = daf(:);
dbf = dbf(:);
y = [daf(1:dmid)' dbf(dmid+1:end)'];


end