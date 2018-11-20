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
