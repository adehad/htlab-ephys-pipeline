%%%% SpikeAnalysis3
%%%% A. Haddad, H. Lin
%%%% INPUT: 
%%%% Combine code segments from spikegl_harvest_HTL.m, SpikeRastor.m

%% Check Required Functions
addpath(genpath('requiredFunctions')) % path to folder with functions
reqFuncStr = '';
if ~exist('readMetafile.m','file'); reqFuncStr=strcat(reqFuncStr,'readMetafile.m, ');           end
if ~exist('get_ch_thr.m','file');   reqFuncStr=strcat(reqFuncStr,'get_ch_thr.m, ');             end
if ~exist('get_events.m','file');   reqFuncStr=strcat(reqFuncStr,'get_events.m (get_ch_thr), ');end
if ~exist('splitconv.m','file');    reqFuncStr=strcat(reqFuncStr,'splitconv.m (get_events), '); end
if ~exist('tconv.m','file');        reqFuncStr=strcat(reqFuncStr,'tconv.m (splitconv), ');      end
if ~exist('wavecull.m','file');     reqFuncStr=strcat(reqFuncStr,'wavecull.m, ');               end
if ~exist('rasterplot.m','file');   reqFuncStr=strcat(reqFuncStr,'rasterplot.m, ');             end
if ~exist('plot_dir_ade.m','file');  reqFuncStr=strcat(reqFuncStr,'plot_dir_ade.m (can replace with plot), ');            end
if ~exist('plotellipse.m','file');  reqFuncStr=strcat(reqFuncStr,'plotellipse.m, ');            end


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
filename.raw=['181017_01.bin']; % Raw data file
filename.kilosortOutput=['181017_01_sorted.mat']; %including structure s and trial_spikes, trial_cluster
% [m, fpath, mfile] = readMetafile;
% m.metafile = mfile;
% m.metapath = fpath;
% NO MetaFile with Daniel's Data, manually entering:
m.StimGL_nloops = 20;
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
m.pdthr = get_ch_thr(filename.raw,m.nChans,m.pdch,5e5,'pd',m.pd_f,3,0,0);
% m.pdthr = 3;

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

%% start to construct data for rastor plot
%%%% INPUT: waves,trial_spikes,pd(above) m.xytraj  m.StimGL_nloops

getxy = 1;
% Get Stimulus Data
% C:\Users\Adehad\Desktop\Box Sync\DF_Signal\Reference\HT spike sorting\Stimuli for 131210
if getxy
    [m.xyname, m.xypath] = uigetfile('*.*','Select stimGL framevar (xytraj) file:');
    VAR=importdata([m.xypath,m.xyname]);
    m.xytraj = VAR.data(:,5:6); %int16(xytraj);
end

% Basic info
D=34;        % 55; Distance of the DF head to the screen
W=125;       % 123; Width of the projection
H=70;        % Height of the projection
xMap=1280/W; % 111; % pix/mm 125mm
yMap=720/H;  % 63; % pix/mm % 80mm
dt=1/360;    % 1/fps    fps is of projector
stimLength=size(m.xytraj,1);
sort_bit=0; %to plot raster and RF for sorted and non-sorted data

m.pd = pd;
% waves=double(waves);            % waves already created above
spikes=double(trial_spikes);    % trial_spikes is location in binary file
spikes = double(spikes - m.pd(1) + 1); %Start aligned the spikes  /m.msec;
pd_diff = double(diff(m.pd));

Repeat_thrs = 0.5; % this is the cutoff threshold for how consistent the spikes are over repeated trials. 



%% generate receptive field plot

%% Step 1: Reconstruct target angular data
elev_offset=300; %600;
xytraj_temp(:,1)=-double(m.xytraj(:,1))+1280; % invert the image back
xytraj_temp(:,2)=-double(m.xytraj(:,2))+720;
TarTraj(:,1)=atan((xytraj_temp(:,1)-640)/(D*xMap))*180/pi;
TarTraj(:,2)=atan((xytraj_temp(:,2)-720+elev_offset)/(D*yMap))*180/pi;
TarTraj_dot=diff(TarTraj,1)/dt;
TarTraj_dot(:,1)=csaps(1:stimLength-1,TarTraj_dot(:,1),0.5,1.5:stimLength-0.5);
TarTraj_dot(:,2)=csaps(1:stimLength-1,TarTraj_dot(:,2),0.5,1.5:stimLength-0.5);
TarTraj_dot=[TarTraj_dot(1,:); TarTraj_dot];

%% Step 2: Align all the data series
m.pd_diff_threshold = 1.6e3; % during each repeat of the stimulus there is a repeated PD event, this threshold finds it
[RepeatSpace, RepeatIndex]=find(pd_diff>m.pd_diff_threshold); % this is in 360Hz unit because diff is actually in 360Hz 
nloops = m.StimGL_nloops;
xytraj = m.xytraj;
% StartFs =StartFrame;
TarTraj_s=TarTraj;
TarTraj_dot_s=TarTraj_dot;
for nl=1:nloops-1
    RepeatTime(nl)=round(double(m.pd(RepeatIndex(nl)+1)-m.pd(1))*360/40000);
    xytraj = [xytraj; zeros((RepeatTime(nl)-size(xytraj,1)),2)];
    xytraj = [xytraj; m.xytraj];
%     StartFs = [StartFs RepeatTime(nl).*ones(1,size(StartFrame,1))+StartFrame];    
    TarTraj_s = [TarTraj_s; zeros((RepeatTime(nl)*-size(TarTraj_s,1)),2)];
    TarTraj_s = [TarTraj_s; TarTraj];
    TarTraj_dot_s = [TarTraj_dot_s; zeros((RepeatTime(nl)-size(TarTraj_dot_s,1)),2)];
    TarTraj_dot_s = [TarTraj_dot_s; TarTraj_dot];    
end

%% Step 3: Compile spike raster

selected_unit = 9; % set to a 0 indexed number, or to 'all' to show all units

if length(selected_unit)>1
    sort_bit = 0;
else 
    sort_bit = 1;
    UnitInd = ismember(trial_spikes, s.(['unit_',num2str(selected_unit)]) );
end

if sort_bit==0
    spikes_360=round(spikes*360/40000);
elseif sort_bit==1
    spikes_360=round(spikes(UnitInd)*360/40000);  
end
RasterMask=zeros(1,stimLength);
ind=find(spikes_360>0 & spikes_360<stimLength);
RasterMask(spikes_360(ind))=ones(1,size(ind,2));

Raster(1,:)=RasterMask; 
SpikeTiming=find(Raster(1,:));
T=RasterMask; % linear raster array
for re=1:nloops-1
    RasterMask=zeros(1,stimLength);
    temp=spikes_360-RepeatTime(re);
    ind=find(temp>0 & temp<stimLength);
    RasterMask(temp(ind))=ones(1,size(ind,2));
%     SpikeTiming=[SpikeTiming temp(ind)];
    SpikeTiming=[SpikeTiming find(RasterMask)];
    Raster(re+1,:)=RasterMask; 
    T=[T RasterMask];
end
TE=find(T==1);
% SpikeTiming=spikes_360(ind);
[nelements,centers]=hist(SpikeTiming,0:36:stimLength);
SpikeHis=csaps(centers,nelements,0.5,1:stimLength); % smooth and upsample


figure(2)
subplot(10,1,1:2); % Position
    plot(TarTraj(:,1),'r'); hold on; plot(TarTraj(:,2),'b');
    xlim([0 stimLength])
    ylim([-60 60])
subplot(10,1,3:4);  % Velocity
    plot(TarTraj_dot(:,1),'r'); hold on; plot(TarTraj_dot(:,2),'b');
    xlim([0 stimLength])
    ylim([-150 150])

if nloops>1
    ntrials = nloops;
else
    ntrials = nloops;
end
subplot(10,1,5:9);  % Spike Raster
    rasterplot(TE,ntrials,stimLength,gca);
%     ylabel(['Trails (',num2str(ntrials),')']) % incorporated to rasterplot
    % rasterplot(TE,nloops-1,stimLength,'plotwidth',0.5);
subplot(10,1,10);   % Spike Hist
    %bar(centers./360,nelements);
    plot((1:stimLength)./360,SpikeHis)
    xlim([0 stimLength./360])
    ylim([0 50])

%% Save post-Rastor
% save([filename.kilosortOutput(1:10),'_postRaster'])

%% Plot Receptive Fields
%%%% INPUT: TarTraj (above), Gaze_rest, head_offset, spikes_180 (spike times in 180Hz projector indices)
%%%% Requires: plotellipse.m
%% set up some variables
head_offset_ang=0;  %this depends on the head marker mounting azimuthal accuracy per animal
Gaze_rest=51;       %resting gaze elevation in degree
spikes_180 = spikes_360';
rejected_spikes_180 = spikes_180(spikes_180 <= 0 );
spikes_180 = spikes_180(spikes_180 > 0 );

disp(['# Negative spikes omitted = ', num2str(length(rejected_spikes_180))]);


% saccade cone parameters:
ab=44*pi/180; %88
bb=56.5*pi/180; %113/2=
zb=[0 75*pi/180];%estimated saccade cone axis 75deg
alphab=0;

% takeoff cone parameters:
ab2=22.5*pi/180; %45
bb2=28.5*pi/180; %57
zb2=[0 80*pi/180];%estimated saccade cone axis 80deg
alphab2=0;
spike_tail = 5; % how many samples to plot leading up to the spike location
response_offset = 5; % how many samples to shift spike location in receptive field to account for animal response time
%% plot spike time stimuli in angular coordinate
figure(3); hold on;
ax_fig3 = gca;
plot(TarTraj(:,1), TarTraj(:,2)+Gaze_rest+head_offset_ang,'.','color',[0.8 0.8 0.8]); 
errorInLoop = [];

% Fix error in KiloSort - can comment out once fixed
spikes_180 = unique(spikes_180);

for ss=1:ntrials
    
    spikes_inTrajTrial = ( spikes_180>( (ss-1)*(size(TarTraj,1)) ) ) & ( spikes_180<( ss*(size(TarTraj,1)) ) ) ;
       
    spikes_inTrajTrial = spikes_180(spikes_inTrajTrial); %
    
    spikes_inTrajTrial = spikes_inTrajTrial - response_offset;
    
    spikes_inTraj.(['trial_',num2str(ss)]) = spikes_inTrajTrial ;
        % Stores what spikes were found for a given repeat of stimuli
        
    
    fprintf('Plotting Trial: %02d, %03d (Spikes) \n', ss, length(spikes_inTrajTrial));
    title(['Plotting Trials 1 to ' , num2str(ss)])
    for kk=1:length(spikes_inTrajTrial)
        temp_startElem = spikes_inTrajTrial(kk);
        temp_tailElem = temp_startElem - spike_tail;
        if temp_tailElem < (ss-1)*(size(TarTraj,1))
            temp_tailElem = (ss-1)*(size(TarTraj,1))+1;
        end
        
%--------% Spike Location
         plot(TarTraj_s(temp_startElem,1), ...
              TarTraj_s(temp_startElem,2)+Gaze_rest+head_offset_ang,'r.');
     
%--------% Spike tail
        xTrajVec = TarTraj_s( temp_tailElem:temp_startElem, 1 );
        yTrajVec = TarTraj_s( temp_tailElem:temp_startElem, 2 )+Gaze_rest+head_offset_ang;
        
        % filter out any movement that is too large
            movementThresh = 12; % in angular degrees
            
            validVectorElems = abs(diff(xTrajVec))<movementThresh & abs(diff(yTrajVec))<movementThresh;

            temp_EndElem = find(validVectorElems==true,1,'last');
            if isempty(temp_EndElem) % empty if there are no valid elements
                temp_EndElem = length(validVectorElems);
            end

            temp_StartElem = find(validVectorElems(1:temp_EndElem)==false,1,'last');
            if isempty(temp_StartElem) % empty if all remaining are valid
                temp_StartElem = 0; % will +1 below
            end

            validVectorElems = (temp_StartElem+1:temp_EndElem+1);
        try
        xTrajVec = xTrajVec(validVectorElems);
        yTrajVec = yTrajVec(validVectorElems);
        catch ME
            warning(ME.message)
        end
%         plot( ax_fig3, xTrajVec, yTrajVec, 'b' ); % Faster but less useful I think
        plot_dir_ade( xTrajVec, ...            
                      yTrajVec, [], ax_fig3);
%           drawnow;   % use the one outside this for loop to be faster  
    end
    drawnow;  % use this if you want to see the spikes per trial live
end

 
plotellipse(zb*180/pi, bb*180/pi, ab*180/pi, alphab, 'c--')
plot(zb(1)*180/pi,zb(2)*180/pi,'c+')

plotellipse(zb2*180/pi, bb2*180/pi, ab2*180/pi, alphab2, 'm--')
plot(zb2(1)*180/pi,zb2(2)*180/pi,'m+')

plot(0,Gaze_rest,'k+'); axis equal
    xlim([-60 60]); ylim([0 120])
    xlabel('azimuth (deg)'); ylabel('elevation (deg)'); 
    title('Target in body oriented global ref')
figure(3)
hold off
%% plot spike time stimuli in angular coordinate [Heat Map]
% Using TarTraj_s - as this includes the repeats
figure(4); hold on;
% scatter(TarTraj(:,1), TarTraj(:,2)+Gaze_rest+head_offset_ang,10,[0.8 0.8 0.8],'filled'); 
% alpha 0.5
numBins = 6; % Number of bins for each dimensions
histogram2( TarTraj_s(spikes_180(spikes_180<=length(TarTraj_s)),1), ...
            TarTraj_s(spikes_180(spikes_180<=length(TarTraj_s)),2)+Gaze_rest+head_offset_ang, ...
            numBins, ...
            'DisplayStyle','tile','ShowEmptyBins','on' )
% Omits spikes that are outside the length of the stimuli

plotellipse(zb*180/pi, bb*180/pi, ab*180/pi, alphab, 'c--')
plot(zb(1)*180/pi,zb(2)*180/pi,'c+')

plotellipse(zb2*180/pi, bb2*180/pi, ab2*180/pi, alphab2, 'm--')
plot(zb2(1)*180/pi,zb2(2)*180/pi,'m+')

plot(0,Gaze_rest,'k+'); axis equal
    xlim([-60 60]); ylim([0 120])
    xlabel('azimuth (deg)'); ylabel('elevation (deg)'); 
    title('Target in body oriented global ref')

%% overlay a head map of spike density
