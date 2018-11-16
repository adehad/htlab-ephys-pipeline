function extractTrialUnitWavesFileName = extractTrialUnitWaves(rawBinary,extractUnitsFileName)
% Copy Paste JOB for now
%%%% SpikeAnalysis3
%%%% A. Haddad, H. Lin
%%%% INPUT: 
%%%% Combine code segments from spikegl_harvest_HTL.m, SpikeRastor.m
extractTrialUnitWavesFileName = '';
%% Check Required Functions
addpath(genpath('requiredFunctions')) % path to folder with functions
reqFuncStr = '';
if ~exist('readMetafile2.m','file'); reqFuncStr=strcat(reqFuncStr,'readMetafile2.m, ');           end
if ~exist('get_ch_thr.m','file');   reqFuncStr=strcat(reqFuncStr,'get_ch_thr.m, ');             end
if ~exist('get_events.m','file');   reqFuncStr=strcat(reqFuncStr,'get_events.m (get_ch_thr), ');end
if ~exist('splitconv.m','file');    reqFuncStr=strcat(reqFuncStr,'splitconv.m (get_events), '); end
if ~exist('tconv.m','file');        reqFuncStr=strcat(reqFuncStr,'tconv.m (splitconv), ');      end
if ~exist('wavecull.m','file');     reqFuncStr=strcat(reqFuncStr,'wavecull.m, ');               end

if ~strcmp(reqFuncStr,'')
error(['The following functions could not be found in the search path:', newline, reqFuncStr])
end
%% Default Settings
% If the variables do not exist in the workspace, or are empty, set defaults
if ~exist('show','var') || isempty(show);               show = 1;      end
if ~exist('fps','var') || isempty(fps);                  fps = 360;    end
if ~exist('spikewin','var') || isempty(spikewin);   spikewin = 1;      end
if ~exist('pdwin','var') || isempty(pdwin);            pdwin = 4;      end
if ~exist('block','var') || isempty(block);            block = 1e6;    end
if ~exist('esize','var') || isempty(esize);            esize = 16;     end
if ~exist('getxy','var') || isempty(getxy);            getxy = 0;      end

%% set up some parameters in structure m
    % Path can also be added to filenames below
filename.raw=['181017_09.bin']; % Raw data file
filename.kilosortOutput=['181017_09_sorted.mat']; %including structure s and trial_spikes, trial_cluster
% [m, fpath, mfile] = readMetafile;
% m.metafile = mfile;
% m.metapath = fpath;
% NO MetaFile with Daniel's Data, manually entering:
m.StimGL_nloops = 5;
m.nChans = 2;
m.sRateHz = 30e3;
% m.fileSizeBytes
    fid = fopen(filename.raw);
    if fid > 0
    fseek(fid,0,'eof');
    m.fileSizeBytes = ftell(fid);
    fseek(fid,0,'bof');
    else
        error('Could not open binary file, filename.raw')
    end



m.fps       = 360; %projector frame rate (RGB)
m.pdch      = m.nChans; %assume pd is last ch
m.ech       = 1:m.nChans-1; % ephys channel(s) is everything except the last
m.dbytes    = 2; % byte size of data - i.e. int16 is 2 bytes
m.msec      = m.sRateHz/1000; % conversion factor from ms time to sample number
m.spikewin  = spikewin*m.msec; % window to crop around spikes
m.pdwin     = pdwin*m.msec; % window to crop around pd events
esize_msec  = esize*m.msec; % in units of sample pts
esize_bytes = esize_msec*m.nChans*m.dbytes; % in units of bytes
dlen        = m.fileSizeBytes/(m.dbytes*m.nChans);
trueblock   = block - 2*esize_msec;     % Block reading code has fseek the moves back by 2*esize_msec

%% extract waves according to spike times to confirm that the sorting is correct
% Load select KiloSort data
load(filename.kilosortOutput) %including structure s and trial_spikes, trial_cluster

% denoising filters
m.el_flen = 400;
m.pd_flen = 200;
m.el_cut = 400;
m.pd_cut = m.fps/2;
m.el_f = fir1(m.el_flen,[m.el_cut 7000]./(m.sRateHz/2)); % spikes
m.pd_f = fir1(m.pd_flen,m.pd_cut/(m.sRateHz/2),'low'); % photodiode


% set event thresholds
% m.thr = get_ch_thr(filename.raw,m.nChans,m.ech,5e5,'electrode',m.el_f,1e3,1,0);
% m.pdthr = get_ch_thr(filename.raw,m.nChans,m.pdch,6e5,'pd',m.pd_f,3,0,0);
m.pdthr = 3;

x = [];
spikes = [];
waves = [];             
incompleteWave = [];
pd = [];
pwaves = [];
slong_srate = 1000; %Hz
slong_drate = round(m.sRateHz/slong_srate);
m.slong_srate = m.sRateHz/slong_drate;
slong_block = trueblock/slong_drate;
slong_block = ceil(trueblock/slong_drate); % Need to do this as dimension mismatch
m.slong = int16(zeros(1,slong_block*floor(dlen/trueblock)));

blocklen_k = 1e6;
fid = fopen(filename.raw,'r');
lastblock = floor(dlen/trueblock) + 1;
% figure(111); hold off; clf
for k = 1:lastblock  % go to the last block
    if k == lastblock
        blocklen_k = inf;
    else
        blocklen_k = block;
    end
    x = fread(fid,[m.nChans blocklen_k],'int16'); % load one block at a time
 
    % get downsampled raw data for long time series examination; cull out
    % block overlaps appropriately. NOTE: leaves out the last piece of time
    % series (that is < block pts in length) !!
    if k < lastblock
        m.slong(1+(k-1)*slong_block:k*slong_block) = decimate(x(m.ech,1:end-2*esize_msec),slong_drate);
    end
    
    % get pd data
    [pdk,pwavetmp] = get_events(x,m.pdch,m.pd_f,m.pdwin,m.msec,m.pdthr,0,0);
    truepd = find((pdk > esize_msec) & (pdk < size(x,2) - esize_msec));
    pwaves = [pwaves; pwavetmp(truepd,:)];
    pd = [pd (pdk(truepd)+(k-1)*(trueblock))];

    
    % get electrode data from KiloSort output files

    % First fix any incomplete spikes with POSITIVE indexes
    if ~isempty(incompleteWave)
        for inc = 1:1:size(incompleteWave,1)
            if incompleteWave(inc,2) > 0
                waves(:, end-incompleteWave(inc,2):end, incompleteWave(inc,1)) = x(m.ech,end-incompleteWave(inc,2):end);
                    % takes the first elements of the newly loaded x and
                    % stores into the last element of the wave (inc, 2)
            end
        end
        incompleteWave = incompleteWave(incompleteWave(:,2)<0, :);
        % remove all incompleteWaves that have been fixed leaving behind
        % only the NEGATIVE indexes
    end
    
    if k < lastblock
        b_Offset = (k-1)*blocklen_k;
        fseek_Offset = (k-1)*(+2*esize_msec); % Set to 0 if you are not using the fseek -2*esize_bytes at the end of the block loop
        b_Offset = b_Offset - fseek_Offset;
    else % Comment this else part out if you do not want to search for spikes in the final block
        blocklen_k = block;
        b_Offset = (k-1)*blocklen_k;
        fseek_Offset = (k-1)*(+2*esize_msec); % Set to 0 if you are not using the fseek -2*esize_bytes at the end of the block loop
        b_Offset = b_Offset - fseek_Offset;
        blocklen_k = size(x,2);
    end
    
    row = []; col =[]; b_spikes =[]; % Reset values every loop
    [row, col, b_spikes]= find(trial_spikes > (b_Offset+(2*esize_msec)) & trial_spikes <= (b_Offset+blocklen_k));    
                % trial_spikes equivalent to pdk (event times)
                % +(2*esize_msec)... is because blocks overlap due to fseek - we don't want to record same spike multiple times
                
    tempWaves = zeros( length(m.ech), 2*m.spikewin+2, length(b_spikes) ); % #channels by #elements per channel by #spikes
        % Store how each spike looks in all the data channels
        
    tempIncompleteWaves = int64(zeros(length(b_spikes), 3)); % [ #spikes (1 if wave is incomplete), (stating index of incomplete wave), (stating how many missing elements) ]
        % Store Info about incomplete waves
        
    for ss=1:length(b_spikes) % go through each spike in the block and assign them into units; extract spike waveform +-1ms around peak
        tempEl = trial_spikes(row(ss),col(ss));
            % Temporary storage of element in trial_spikes - i.e. location of spike
        tempIdx = [int64(tempEl)-b_Offset-m.spikewin:int64(tempEl)-b_Offset+m.spikewin+1];
            % Temporary indexes storing where in the loaded data x, to take
            % waveform data , +1 so even number 
            
        % Handle when Index is out of bounds of loaded data
        if tempIdx(1)<1     % negative index
            tempWave = zeros(length(m.ech), length(tempIdx));
            tempWave(:, end-length(tempIdx(tempIdx>0))+1:end) = x(m.ech,tempIdx(tempIdx>0)); 
             % Stores only elements in range into tempWave
             % x data is # channels by # elements (per chan) - hence need to transpose to match tempWave
            tempWaves(:,:,ss) = tempWave;
%             waves = [waves; tempWave];  
            tempIncompleteWaves(ss,:) = [1, size(waves,3)+ss, int64(tempIdx(1))];
%             incompleteWave = [incompleteWave; size(waves,1), int64(tempIdx(1))];
                % store which wave is incomplete and how many missing
                % elements
        elseif tempIdx(end)>blocklen_k  % index greater than block
            tempWave = zeros(length(m.ech), length(tempIdx));
            tempWave(:, 1:length(tempIdx(tempIdx<=blocklen_k))) = x(m.ech,tempIdx(tempIdx<=blocklen_k));
             % Stores only elements in range into tempWave
             % x data is # channels by # elements (per chan) - hence need to transpose to match tempWave
             tempWaves(:,:,ss) = tempWave;
%             waves = [waves; tempWave];  
            tempIncompleteWaves(ss,:) = [1, size(waves,3)+ss, int64(tempIdx(end))-int64(blocklen_k)];
%             incompleteWave = [incompleteWave; size(waves,1), int64(tempIdx(end))-int64(blocklen_k)];
                % store which wave is incomplete and how many missing
                % elements
        else
            tempWave = zeros(length(m.ech), length(tempIdx));
            tempWave(:, :) = x(m.ech,tempIdx);
            tempWaves(:,:,ss) = tempWave;
%             waves = [waves; x(m.ech,tempIdx)];
        end

    end
    
    %%%%%%%%%%%%%%%%%%%%%%
    waves = cat(3,waves, tempWaves); 
        % Join arrays at the 3rd dimension
    tempIncompleteWaves = tempIncompleteWaves( tempIncompleteWaves(:,1)~=0, :);
        % Keep only waves that are actually incomplete, column 1 =1 if wave incomplete
    incompleteWave = [incompleteWave; tempIncompleteWaves(:,2:3)];
        % Columns 2 and 3 contain the useful information
    
    
        % First fix any incomplete spikes with NEGATIVE indexes
    if ~isempty(incompleteWave)
        for inc = 1:size(incompleteWave,1)
            if incompleteWave(inc,2) < 0
                waves(:, 1:abs(incompleteWave(inc,2))+1, incompleteWave(inc,1)) = negativeX(m.ech,1:abs(incompleteWave(inc,2))+1);
                    % +1 to fix a problem if value is 0
                    % negativeX is already flipped so just takes the same
                    % number of elements as specified by (inc,2)
            end
        end
        incompleteWave = incompleteWave(incompleteWave(:,2)>0, :);
            % remove all incompleteWaves that have been fixed, leaving
            % behind only the POSITIVE indexes
    end
    negativeX = fliplr(x(m.ech,end-m.spikewin+1:end)); % storage buffer to fix incomplete spikes

    fprintf('Data block %02d of %d; %d spikes found, %d pd events \n',k,floor(dlen/block)+1, ...
                                    size(tempWaves,3), length(pdk))
    fseek(fid,-2*esize_bytes,'cof');
%     figure(111); plot(squeeze(waves(:,:,(size(waves,3)-size(tempWaves,3)+1):50:end))); title(['Data Blocks 1 to ', num2str(k)]); hold on; drawnow; 
end
fclose(fid);
disp('-----------------------------------------------------------')

% Convert waves structure
% Input:  rows X columns X pages = channels X spike elements X spike e.g
    % (1,:,1) - channel 1, all elements of the first spike
% Output: rows X columns X pages = spikes   X spike elements X channels 
    % (1,:,1) - first spike, all elements of channel 1
waves = permute(waves, [3 2 1]); formatChanged = true;
    % 1,2,3 refer to dimension, so we have switched dimensions 1 and 3

 

%% Add waveforms to the struct s
    % First clear out any existing waves
origFieldNames = fieldnames(s); validFieldNames = [];
validClusterNumbers = [];
for ii=1:length(origFieldNames) 
    if sum(ismember(origFieldNames{ii},'unit_'))>=5
        validFieldNames{end+1} = origFieldNames{ii};
        validClusterNumbers{end+1} = str2num(origFieldNames{ii}(6:end));
    elseif sum(ismember(origFieldNames{ii},'waves_'))>=6
        s = rmfield(s,origFieldNames{ii});
    end
end
    % Now add the waves to the struct
for ii=1:length(validFieldNames)
    unitWaves = ismember(trial_spikes, s.([ 'unit_', num2str(validClusterNumbers{ii}) ]) );
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    if ~formatChanged
        unitWaves = unitWaves(1:size(waves,3)); % For channels X spike elements X spike format
    else
        unitWaves = unitWaves(1:size(waves,1));  % For spikes   X spike elements X channels 
    end
        % Because the code ignores the last chunk of data smaller than the
        % block size, the actual number of spikes found by KiloSort
        % (trial_spikes) is larger than those we process (waves, 3rd
        % dimension)
        % This has now been corrected in the block reading in the else loop that involves: blocklen_k = size(x,2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~formatChanged
        s.([ 'waves_', num2str(validClusterNumbers{ii}) ]) = waves(:,:,unitWaves); % For channels X spike elements X spike format
    else
        s.([ 'waves_', num2str(validClusterNumbers{ii}) ]) = waves(unitWaves,:,:); % For spikes   X spike elements X channels 
    end
end

%% Save pre-Rastor
% save([filename.kilosortOutput(1:10),'_preRaster'])



end