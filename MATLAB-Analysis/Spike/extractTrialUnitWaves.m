function [m,s] = extractTrialUnitWaves(rawBinary, extractUnitsFileName, m, outputFileName )
%%%% A. Haddad, H. Lin
%%%% INPUT: 
%%%%     % Path can also be added to filenames of rawBinary / extractUnitsFileName
% if outputFileName = [], it will not save this function's data as a mat
%% Check Required Functions
% Appended to the end of this fuction
%       [None]
%% Default Settings
spikewin = 1;        % size in ms to crop around spikes  
block = 1e6;         % block size for loading data
esize = 16;          % ?

%% set up some parameters in structure m
esize_msec  = esize*m.msec; % in units of sample pts
esize_bytes = esize_msec*m.nChans*m.dbytes; % in units of bytes
dlen        = m.fileSizeBytes/(m.dbytes*m.nChans);
trueblock   = block - 2*esize_msec;     % Block reading code has fseek the moves back by 2*esize_msec

m.spikewin  = spikewin*m.msec; % window to crop around spikes


%% extract waves according to spike times to confirm that the sorting is correct
% Load select KiloSort data
load(extractUnitsFileName) %including structure s and trial_spikes, trial_cluster

% denoising filters
m.el_flen = 400;
m.el_cut = 400;
m.el_f = fir1(m.el_flen,[m.el_cut 7000]./(m.sRateHz/2)); % spikes


x = [];
waves = [];             

blocklen_k = 1e6;
fid = fopen(rawBinary,'r');
if fid > 0
    fseek(fid,0,'eof');
    m.fileSizeBytes = ftell(fid);
    fseek(fid,0,'bof');
else
    fclose(fid);
    error(['Could not open binary file: ', rawBinary])
end

lastblock = floor(dlen/trueblock) + 1;
spikesFound = [];
% figure(111); hold off; clf
for k = 1:lastblock  % go to the last block
    
    if k == lastblock % the final block is smaller than blocklen_k, inf, will load everything it can
        blocklen_k = inf;
    else
        blocklen_k = block;
    end
    
    x = fread(fid,[m.nChans blocklen_k],'int16'); % load one block at a time
    
    for cc=1:length(m.ech)  % Filder Ephys Channels
        x(cc,:) = splitconv(x(m.ech(cc),:),m.el_f);
    end

    % get electrode data from sorted output files
    
    % Establish b_Offset - the starting element of the block - so we can search trial_spikes properly 
    if k < lastblock
        b_Offset = (k-1)*blocklen_k;
        fseek_Offset = (k-1)*(+2*esize_msec); 
        b_Offset = b_Offset - fseek_Offset;
    else % Comment this else part out if you do not want to search for spikes in the final block
        blocklen_k = block;
        b_Offset = (k-1)*blocklen_k;
        fseek_Offset = (k-1)*(+2*esize_msec); 
        b_Offset = b_Offset - fseek_Offset;
        blocklen_k = size(x,2);
    end
    
    
    row = []; col =[]; b_spikes =[]; % Reset values every loop
    % Ensure we get the most number of extractable spikes - by making sure
    % spike index has enough space to extract waveform
    if k ==1
        [row, col, b_spikes]= find( trial_spikes > (b_Offset+m.spikewin+1) & trial_spikes <= (b_Offset+blocklen_k-esize_msec) );
            % First block should only look after m.spikewin to prevent
            % negative spikes
    elseif k==lastblock
        [row, col, b_spikes]= find( trial_spikes > (b_Offset+esize_msec) & trial_spikes <= (b_Offset+blocklen_k-m.spikewin-1) );    
            % Last block should look until m.spikewin before the end to
            % prevent positive spikes
    else
        [row, col, b_spikes]= find( trial_spikes > (b_Offset+esize_msec) & trial_spikes <= (b_Offset+blocklen_k-esize_msec) );    
                % trial_spikes equivalent to pdk (event times)
                % +/-(esize_msec)... is because blocks overlap due to fseek - we don't want to record same spike multiple times
    end
    
    tempWaves = zeros( length(m.ech), 2*m.spikewin+2, length(b_spikes) ); % #channels by #elements per channel by #spikes
        % Store how each spike looks in all the data channels
        
    spikesFound = [spikesFound;row]; % keep track of what spikes we have been able to extract
    for ss=1:length(b_spikes) % go  through each spike in the block and assign them into units; extract spike waveform +-1ms around peak
        tempEl = trial_spikes(row(ss),col(ss));
            % Temporary storage of element in trial_spikes - i.e. location of spike
        tempIdx = [ (int64(tempEl)-b_Offset-m.spikewin):(int64(tempEl)-b_Offset+m.spikewin+1) ];
            % Temporary indexes storing where in the loaded data x, to take
            % waveform data , +1 so even number 
            
        % Handle when Index is out of bounds of loaded data - but should
        % never reach here - due to if loop for find( trial_spikes .. ) above
        if tempIdx(1)<1     % negative index
            % do nothing - corrected for by fseek() shift
            disp('a')
        elseif tempIdx(end)>blocklen_k  % index greater than block
            % do nothing - corrected for by fseek() shift
            disp('b')
        else
            tempWave = zeros(length(m.ech), length(tempIdx));
            tempWave(:, :) = x(m.ech,tempIdx);
            tempWaves(:,:,ss) = tempWave;
        end

    end
    
    %%%%%%%%%%%%%%%%%%%%%%
    waves = cat(3,waves, tempWaves); 
        % Join arrays at the 3rd dimension

    fprintf('Data block %02d of %d; %d spikes found \n' ...
            ,k,floor(dlen/block)+1, size(tempWaves,3))
                                
    % Shift back - creates overlapping blocks
    fseek(fid,-2*esize_bytes,'cof');
        
%     figure(111); plot(squeeze(waves(:,:,(size(waves,3)-size(tempWaves,3)+1):50:end))); title(['Data Blocks 1 to ', num2str(k)]); hold on; drawnow; 
end
fclose(fid);
disp('-----------------------------------------------------------')

% There will be some spikes on either end of the dataset where we cannot
% extract waveforms - the spikewin will overflow the original binary
disp([num2str(length(trial_spikes)-size(waves,3)),' spikes are unextractable']);

% Add zeros to account for this
waves = cat(3, zeros(length(m.ech), length(tempIdx),spikesFound(1)-1), waves); % To beginning of waves
waves = cat(3, waves, zeros(length(m.ech), length(tempIdx),size(trial_spikes,1)-spikesFound(end))); % To end of waves


formatChanged = false;
% Convert waves structure
% Input:  rows X columns X pages = channels X spike elements X spike e.g
    % (1,:,1) - channel 1, all elements of the first spike
% Output: rows X columns X pages = spikes   X spike elements X channels 
    % (1,:,1) - first spike, all elements of channel 1
waves = permute(waves, [3 2 1]); formatChanged = true;
    % 1,2,3 refer to dimension, so we have switched dimensions 1 and 3


%% Add waveforms to the struct s
fprintf('Saving waveforms  ... ');
    % First clear out any existing waves
origFieldNames = fieldnames(s); validFieldNames = [];
validClusterNumbers = [];
for ii=1:length(origFieldNames) 
    if contains(origFieldNames{ii},'unit_')
        validFieldNames{end+1} = origFieldNames{ii};
        validClusterNumbers{end+1} = origFieldNames{ii}(6:end); %str2num(origFieldNames{ii}(6:end));
    elseif contains(origFieldNames{ii},'waves_')
        s = rmfield(s,origFieldNames{ii});
    end
end
    % Now add the waves to the struct
for ii=1:length(validFieldNames)
    unitWaves = ismember(trial_spikes, s.([ 'unit_', validClusterNumbers{ii} ]) );
    if formatChanged
        s.([ 'waves_', num2str(validClusterNumbers{ii},'%02i') ]) = waves(unitWaves,:,:); % For spikes   X spike elements X channels
    else
        s.([ 'waves_', num2str(validClusterNumbers{ii},'%02i') ]) = waves(:,:,unitWaves); % For channels X spike elements X spike format
    end
end
fprintf('done! \n');

%% Save to .mat File

if ~isempty(outputFileName) % only save if filename is not empty
    fprintf('Saving function output ... ');
    % clear up variables before saving to keep it neat
    clear fid ans trial_clusters tempEl tempWave x tempIdx tempWaves formatChanged unitWaves validFieldNames validClusterNumbers origFieldNames row col b_spikes blocklen_k ii ss k cc fseek_Offset b_Offset lastblock
    save(outputFileName)
    fprintf('done! \n');
end
end