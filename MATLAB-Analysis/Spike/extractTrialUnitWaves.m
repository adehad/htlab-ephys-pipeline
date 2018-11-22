function [m,s] = extractTrialUnitWaves(rawBinary, extractUnitsFileName, m, spike_screen_bit, outputFileName )
%%%% A. Haddad, H. Lin
%%%% INPUT: 
%%%%     % Path can also be added to filenames of rawBinary / extractUnitsFileName
% spike_screen_bit  - 1: if you want to do secondary template matching
% if outputFileName = [], it will not save this function's data as a mat
%% Check Required Functions
% Appended to the end of this fuction
%       [None]
%% Default Settings
spikewin = 1;        % size in ms to crop around spikes  
block = 1e6;         % block size for loading data
esize = 16;          % ?

% spike_screen_bit = 0;

%% set up some parameters in structure m

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

esize_msec  = esize*m.msec; % in units of sample pts
esize_bytes = esize_msec*m.nChans*m.dbytes; % in units of bytes
dlen        = m.fileSizeBytes/(m.dbytes*m.nChans);
trueblock   = block - 2*esize_msec;     % Block reading code has fseek the moves back by 2*esize_msec

m.spikewin  = spikewin*m.msec; % window to crop around spikes


%% Extract waves according to sorted output's spike times
% Load select KiloSort data
    load(extractUnitsFileName) %including structure s and trial_spikes, trial_cluster
    % First clear out any existing waves
    origFieldNames = fieldnames(s); validFieldNames = []; % names in the struct s
    validClusterNumbers = [];
    for ii=1:length(origFieldNames) 
        if contains(origFieldNames{ii},'unit_')                     % store fieldnames that start with unit_
            validFieldNames{end+1} = origFieldNames{ii};
            validClusterNumbers{end+1} = origFieldNames{ii}(6:end); % stores cluster number
        elseif contains(origFieldNames{ii},'waves_')                % fieldnames that start with waves_ are deleted
            s = rmfield(s,origFieldNames{ii});
        end
    end
    % Redefine trial_spikes (in s) to spike_times - only include selected units spike_times
    spike_times = [];
    for ii=1:length(validFieldNames)
        spike_times = cat(1, spike_times, s.([ 'unit_', validClusterNumbers{ii} ]));
    end
    spike_times = sort(spike_times); % sorts in ascending order - i.e. increasing spike times
    
% denoising filters
m.el_flen = 400;
m.el_cut = 400;
m.el_f = fir1(m.el_flen,[m.el_cut 7000]./(m.sRateHz/2)); % spikes


x = [];
waves = [];             

blocklen_k = 1e6;

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
    
    for cc=m.ech(1):m.ech(end)  % Filter Ephys Channels
        x(cc,:) = splitconv(x(cc,:),m.el_f);
    end

    % get electrode data from sorted output files
    
    % Establish b_Offset - the starting element of the block - so we can search spike_times properly 
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
        [row, col, b_spikes]= find( spike_times > (b_Offset+m.spikewin+1) & spike_times <= (b_Offset+blocklen_k-esize_msec) );
            % First block should only look after m.spikewin to prevent
            % negative spikes
    elseif k==lastblock
        [row, col, b_spikes]= find( spike_times > (b_Offset+esize_msec) & spike_times <= (b_Offset+blocklen_k-m.spikewin-1) );    
            % Last block should look until m.spikewin before the end to
            % prevent positive spikes
    else
        [row, col, b_spikes]= find( spike_times > (b_Offset+esize_msec) & spike_times <= (b_Offset+blocklen_k-esize_msec) );    
                % spike_times equivalent to pdk (event times)
                % +/-(esize_msec)... is because blocks overlap due to fseek - we don't want to record same spike multiple times
    end
    
    tempWaves = zeros( length(m.ech), 2*m.spikewin+2, length(b_spikes) ); % #channels by #elements per channel by #spikes
        % Store how each spike looks in all the data channels
        
    spikesFound = [spikesFound;row]; % keep track of what spikes we have been able to extract
    for ss=1:length(b_spikes) % go  through each spike in the block and assign them into units; extract spike waveform +-1ms around peak
        tempEl = spike_times(row(ss),col(ss));
            % Temporary storage of element in spike_times - i.e. location of spike
        tempIdx = [ (int64(tempEl)-b_Offset-m.spikewin-1):(int64(tempEl)-b_Offset+m.spikewin) ];
            % Temporary indexes storing where in the loaded data x, to take
            % waveform data , +1 so even number 
            
        % Handle when Index is out of bounds of loaded data - but should
        % never reach here - due to if loop for find( spike_times .. ) above
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
% extract waveforms - the spikewin will overflow out of the bounds of the original binary
disp([num2str(length(spike_times)-size(waves,3)),' spikes are unextractable']);

% Add zeros to account for this
waves = cat(3, zeros(length(m.ech), length(tempIdx),spikesFound(1)-1), waves); % To beginning of waves
waves = cat(3, waves, zeros(length(m.ech), length(tempIdx),size(spike_times,1)-spikesFound(end))); % To end of waves


%% Convert waves structure
% Input:  rows X columns X pages = channels X spike elements X spike e.g
    % (1,:,1) - channel 1, all elements of the first spike
% Output: rows X columns X pages = spikes   X spike elements X channels 
    % (1,:,1) - first spike, all elements of channel 1
waves = permute(waves, [3 2 1]);
    % 1,2,3 refer to dimension, so we have switched dimensions 1 and 3


%% Add waveforms to the struct s
    % Add the waves to the struct - can be overridden by secondary spike sorting below
for ii=1:length(validFieldNames)
    unitWaves = ismember(spike_times, s.([ 'unit_', validClusterNumbers{ii} ]) );

    s.([ 'waves_', num2str(validClusterNumbers{ii},'%02i') ]) = waves(unitWaves,:,:); % For spikes   X spike elements X channels
    s.waves{validClusterNumbers{ii}} = waves(unitWaves,:,:);
end

%% Perform Secondary Spike Sorting
if spike_screen_bit == 1
    fprintf('Starting secondary spike sorting ... \n');
    s_tempSpikes=[];    s_tempWaves=[];    % temporary variables storing s.units_ and s.waves_ - these will overwride those currently in s
    % Go through each spike AGAIN to do secondary sorting
    for ii=1:length(validFieldNames)
        s_tempSpikes=[];
        tempWaves = s.([ 'waves_', validClusterNumbers{ii} ]); % Compute template from these waveforms
        for ss=1:size(tempWaves,1) 
            trace=tempWaves(ss,:,:); % Select one spike waveform at a time
            trace = permute(trace, [3 2 1]); % permute to row=electrode sites, column=timepoint, page=spikes
            [matches_ind, match_extract] = TemplateMatch3(trace, tempWaves, m.sRateHz);
            
            s_tempSpikes=cat(1,s_tempSpikes, ~isempty(match_extract)); % Create a logical array of whether the template is a good match or not
        end
        s.([ 'unit2_', validClusterNumbers{ii} ]) = boolean(s_tempSpikes); 
        s.unit2{validClusterNumbers{ii}} = boolean(s_tempSpikes); 
            % Only units that match the template are stored as true
        fprintf(['Finished Unit: ', validClusterNumbers{ii}, '\n']);
    end            
end 
fprintf('... done! \n');
    
%% Save to .mat File

if ~isempty(outputFileName) % only save if filename is not empty
    fprintf('Saving function output ... ');
    % clear up variables before saving to keep it neat
    clear fid ans trial_clusters tempEl tempWave x tempIdx tempWaves formatChanged unitWaves validFieldNames validClusterNumbers origFieldNames row col b_spikes blocklen_k ii ss k cc fseek_Offset b_Offset lastblock
    save(outputFileName)
    fprintf('done! \n');
end
end




%%%% TemplateMatch3
%%%% Huai-Ti Lin [Nov 2018]
%%%% This script scan through the input trace for template waves
%%%% INPUT: temp_waves, trace
%%%% OUTPUT: matches2, match_wave2 [index from the trace that matches]
%%%% trace must be longer than the template
%%%% version 2 expects multi-channle trace and template

function [matches2, match_wave2] = TemplateMatch3(trace,temp_waves,sRate,task_type,match_thd2,plot_bit)
if ~exist('match_thd2') || isempty(match_thd2);     match_thd2 = [0.2 0.4]; end % Set threshold for the peak match default  
if ~exist('plot_bit')   || isempty(plot_bit);       plot_bit = 0;           end % Plotting
if ~exist('task_type')  || isempty(task_type);      task_type = 0;          end % 1: for recentering and 0: for no recentering.
if ~exist('sRate')      || isempty(sRate);          sRate = 40e3;           end % Defalt samping rate for NI-DAQ 40k; OpenEphys 30k

t_crop_bit = 1; % 1: crop the template around peak to reduce the matching length
spikeWidth = 0.5; % ms estimated width of an spike - will crop the template to this size around the peak location
%% acquire template info
n_ch=size(trace,1); %check how many channels
if size(temp_waves,1)==1
    template=temp_waves; %if only the template waveform is given
    temp_p2p=max(template,[],2)-min(template,[],2);
    temp_Std=0.2*temp_p2p.*ones(1,size(template,2),n_ch); % assume 20% of peak-to-peak for each channel
else %if a pack of example waveforms were given  
    template=mean(temp_waves,1); %establish template from template waves    
    temp_p2p=max(template,[],2)-min(template,[],2);
    temp_Std=std(temp_waves,0,1); %find the STD of the averaged wave
end

% reformate template into 2D matrix for analysis purpose
template = permute(template, [3 2 1]); template = template(:,:,1);% rows: channels, columns: channel elements
temp_p2p = permute(temp_p2p, [3 2 1]); temp_p2p = temp_p2p(:,:,1);% channels: peak-peak for each channel
temp_Std = permute(temp_Std, [3 2 1]); temp_Std = temp_Std(:,:,1);% rows: channels, columns: std for each channel's elements
[p2p msc]=sort(temp_p2p,'descend'); % p2p actual value of p2p, msc is the original elment in temp_p2p that had that p2p value
% figure; plot(template(:,:)'); pause % debugging to see original template
% crop template
if t_crop_bit ==1
    crop_half = round((spikeWidth/(1000/sRate))/2); %samples to get before and after the peak 
    peakLoc = size(template,2)/2; % by default the extract wave's peak is centered about the peak    
    template = template(:,  peakLoc-crop_half-1:peakLoc+crop_half );    % crop template to size - negative peak still at middle of this array !
%     figure; plot(template(:,:)'); pause % debugging to see croppped template
end

%% match parameters
r1=1; %ref_err=sum(abs((template+(temp_p2p/2))-template));
r2=1; %modulation for match_thd2
match_thd=sum(temp_Std.^2,2)/r2; %set threshold for the deviation match
deviation=[];   
%% picking out the preliminary matches
for tt=1:size(trace,2)-size(template,2)+1; % slide the template along the entire signal trace to find matches
    temp_diff=trace(:,tt:tt+size(template,2)-1)-template;
    deviation(:,tt)=sum(temp_diff.^2,2)/r1; %sum(abs(temp_diff))/r1;    
end

temp_matches ={};
matches=[];
for cc=1:n_ch
    temp_matches{cc}=find(deviation(cc,:)<match_thd(cc));
end
matches=intersect(temp_matches{msc(1)},temp_matches{msc(2)}); %we only require the two most significant channel match
if ~isempty(matches)
    [true_ind]=[1 find(diff(matches)>100)+1];
    matches=matches(true_ind);

    temp_peak=find(template(msc(1),:)==min(template(msc(1),:)));
    preL=temp_peak-1;
    postL=size(template,2)-temp_peak;
    %% refining the matches
    matches2=[];
    match_wave2=[];
    for mm=1:size(matches,2) %go through each spike
        match_wave=trace(:,matches(mm):matches(mm)+size(template,2)-1);
        match_p2p=max(match_wave,[],2)-min(match_wave,[],2);
        %check if the peak matching is good for at least half of the channels
        if sum(abs(min(template,[],2)-min(match_wave,[],2))./temp_p2p <match_thd2(1))>=0.5*n_ch && sum(abs(match_p2p-temp_p2p)./temp_p2p <match_thd2(2))>=0.5*n_ch %check the spike peak     
           if task_type==1
            peak_ind=find(match_wave(msc(1),:)==min(match_wave(msc(1),:)),1); %peak align to the most significant channel
           else
            peak_ind=preL; %   
           end
           matches2=[matches2 matches(mm)+peak_ind]; %compiling the refined matches with wave first peak as the index
           match_wave2 = cat(3,match_wave2, trace(:,matches(mm)+peak_ind-preL:matches(mm)+peak_ind+postL));
        end
    end
    
    if ~isempty(matches2)
        [true_ind]=[1 find(diff(matches2)>100)+1];
        matches2=matches2(true_ind);
        match_wave2 = match_wave2(:,:,true_ind);
        match_wave2 = permute(match_wave2, [3 2 1]);
    end
    
    %% output
    if plot_bit==1
    figure
    subplot(2,2,1); plot(match_wave2(:,:,1)')
    hold on; plot(template(1,:),'k', 'linewidth',2); hold on; ylim([-5000 5000])
    subplot(2,2,2); plot(match_wave2(:,:,2)')
    hold on; plot(template(2,:),'k', 'linewidth',2); hold on; ylim([-5000 5000])
    subplot(2,2,3); plot(match_wave2(:,:,3)')
    hold on; plot(template(3,:),'k', 'linewidth',2); hold on; ylim([-5000 5000])    
    subplot(2,2,4); plot(match_wave2(:,:,4)')
    hold on; plot(template(4,:),'k', 'linewidth',2); hold on; ylim([-5000 5000])
    % pause
    end
else
    matches2=[];
    match_wave2=[];
end


end


%%%% Debug Plot
%%%% AH [Nov 2018]
%%%% Given a list of numbers corresponding to units, will plot before and
%%%% after
%%%% INPUT: temp_waves, trace
%%%% OUTPUT: matches2, match_wave2 [index from the trace that matches]
%%%% trace must be longer than the template
%%%% version 2 expects multi-channle trace and template