function m = extractTrialADC_PD(rawBinary, m, outputFileName )
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
    [pdk,pwavetmp] = get_events(x,m.pdch,m.pd_f,m.pdwin,m.msec,m.pdthr,0,0);
    truepd = find((pdk > esize_msec) & (pdk < size(x,2) - esize_msec));
    pwaves = [pwaves; pwavetmp(truepd,:)];
    pd = [pd (pdk(truepd)+(k-1)*(trueblock))];
    
    %truepdr = find((pdrk > esize_msec) & (pdrk < size(x,2) - esize_msec));
    %pdRise = [pdRise (pdrk(truepdr)+(k-1)*(trueblock))];
    %truepdf = find((pdfk > esize_msec) & (pdfk < size(x,2) - esize_msec));
    %pdFall = [pdFall (pdfk(truepdf)+(k-1)*(trueblock))];
    

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
%m.pdRise = pdRise;
%m.pdFall = pdFall;

%% Save to .mat File

% clear up variables before saving to keep it neat
if ~isempty(outputFileName)
    fprintf('Saving function output ... ');
    clear fid ans pwavetmp x blocklen_k k cc fseek_Offset b_Offset lastblock
    save(outputFileName)
    fprintf('done! \n');
end
end