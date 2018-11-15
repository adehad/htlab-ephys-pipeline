%%%% RF_Analysis - RANDOM WALK
%%%% A. Haddad, H. Lin
%%%% INPUT: 
%%%% Combine code segments from spikegl_harvest_HTL.m, SpikeRastor.m

%% Load Post Raster
load('181017_09__postRaster.mat')
addpath(genpath('requiredFunctions')) % path to folder with functions
%% Plot spike raster
% use selected_unit to identify visual neurons

selected_unit = 3; % set to a 0 indexed number, or to 'all' to show all units

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

%%
response_offset = 5; % how many samples to shift spike location in receptive field to account for animal response time
                     % e.g. if 5, shifts location 5 samples backwards
                     % set to 0 for no offset
figure(6)

spikes_180 = spikes_360';
rejected_spikes_180 = spikes_180(spikes_180 <= 0 );
spikes_180 = spikes_180(spikes_180 > 0 );

disp(['# Negative spikes omitted = ', num2str(length(rejected_spikes_180))]);

spikes_180 = unique(spikes_180);

validSpikes = spikes_180(spikes_180<length(xytraj)); % Find indexes of spikes that occur within xy

xyTrajTrue(:,1) = xytraj(:,2); % x is second row
xyTrajTrue(:,2) = xytraj(:,1); % y is first row

degreesTraj = atan2( diff(xyTrajTrue(response_offset+1:end,2),1), ...
                     diff(xyTrajTrue(response_offset+1:end,1),1) ); %dy/dx
degreesTraj = (degreesTraj)*180/pi;

degreesTraj( find(diff(xyTrajTrue(response_offset+1:end,2),1)== 0 & ...
                  diff(xyTrajTrue(response_offset+1:end,1),1)== 0) ,1 ) = NaN;
    % when there is no change in x or y, we set it to not a number

degreeBin = 45;
    
% ------>   0 degrees
leftToRight = (degreesTraj>0-degreeBin & degreesTraj<0+degreeBin);
    validLRspikes = leftToRight(validSpikes); % Find if these spikes occur during this motion
    xyLtoR = xyTrajTrue(leftToRight,:);    % Find xy traj of this motion
% <------   -180 OR +180 degrees
rightToLeft = (degreesTraj>+180-degreeBin & degreesTraj<=180 | degreesTraj<-180+degreeBin & degreesTraj>=-180);
    validRLspikes = rightToLeft(validSpikes); % Find if these spikes occur during this motion
    xyRtoL = xyTrajTrue(rightToLeft,:);    % Find xy traj of this motion
% /\        +90 degrees
downToUp = (degreesTraj>90-degreeBin & degreesTraj<90+degreeBin);
    validDUspikes = downToUp(validSpikes); % Find if these spikes occur during this motion
    xyDtoU = xyTrajTrue(downToUp,:);    % Find xy traj of this motion
% \/        -90 degrees
upToDown = (degreesTraj>-90-degreeBin & degreesTraj<-90+degreeBin);
    validUDspikes = upToDown(validSpikes); % Find if these spikes occur during this motion
    xyUtoD = xyTrajTrue(upToDown,:);    % Find xy traj of this motion



subplot(2,2,1)
    scatter(xyTrajTrue(leftToRight,1), xyTrajTrue(leftToRight,2), 10,[0.8 0.8 0.8],'filled')
    hold on
    scatter(xyLtoR(validLRspikes,1), xyLtoR(validLRspikes,2), 10,'r','filled')
    xlim([ min(xyTrajTrue(:,1)), max(xyTrajTrue(:,1)) ]); ylim([ min(xyTrajTrue(:,2)), max(xyTrajTrue(:,2)) ]);
    hold off
    title('Left to Right')
subplot(2,2,3)
    scatter(xyTrajTrue(rightToLeft,1), xyTrajTrue(rightToLeft,2), 10,[0.8 0.8 0.8],'filled')
    hold on
    scatter(xyRtoL(validRLspikes,1), xyRtoL(validRLspikes,2), 10,'r','filled')
    xlim([ min(xyTrajTrue(:,1)), max(xyTrajTrue(:,1)) ]); ylim([ min(xyTrajTrue(:,2)), max(xyTrajTrue(:,2)) ]);
    hold off
    title('Right to Left')
subplot(2,2,2)
    scatter(xyTrajTrue(downToUp,1), xyTrajTrue(downToUp,2), 10,[0.8 0.8 0.8],'filled')
    hold on
    scatter(xyDtoU(validDUspikes,1), xyDtoU(validDUspikes,2), 10,'r','filled')
    xlim([ min(xyTrajTrue(:,1)), max(xyTrajTrue(:,1)) ]); ylim([ min(xyTrajTrue(:,2)), max(xyTrajTrue(:,2)) ]);
    hold off
    title('Down to Up')
subplot(2,2,4)
    scatter(xyTrajTrue(upToDown,1), xyTrajTrue(upToDown,2), 10,[0.8 0.8 0.8],'filled')
    hold on
    scatter(xyUtoD(validUDspikes,1), xyUtoD(validUDspikes,2), 10,'r','filled')
    xlim([ min(xyTrajTrue(:,1)), max(xyTrajTrue(:,1)) ]); ylim([ min(xyTrajTrue(:,2)), max(xyTrajTrue(:,2)) ]);
    hold off
    title('Up to Down')
    

%% Heat Maps 
figure(7)
% need to make normalisations possible - probably with anonymous functions
normaliseHist = 'traj'; % 'traj' - normalise according to number of trajectories in the bin used for spikes

colormap(hot) % number inside indicates number of colours to use
numBins = round([ max(xyTrajTrue(:,1))/50, max(xyTrajTrue(:,2))/60]); % Number of bins for each dimensions
                % essentially [480/N, 640/N]        for a 480x640 projection
gaussianSigma = 0; % standard deviation for gaussian for smoothing
                      % set to 0 for no change
                
% imagesc() has all inputs transposed: columns go horizontally - hence are x co-ordinates
               
                % manually added elements are to force histogram to fill plot area
subplot(2,2,1)
[N,c] = hist3([[xyLtoR(validLRspikes,1);0;max(xyTrajTrue(:,1))], ...
              [xyLtoR(validLRspikes,2);0;max(xyTrajTrue(:,2))]], ...
              numBins); N(1,1)=N(1,1)-1; N(end,end) = N(end,end)-1; % subtract elements we added
[N2,c2] = hist3([[xyTrajTrue(leftToRight,1);0;max(xyTrajTrue(:,1))], ...   % Trajectoriess Histogram - i.e. how often is the trajectory in the same bin as we used for spikes
                 [xyTrajTrue(leftToRight,2);0;max(xyTrajTrue(:,2))]], ...
                 numBins); N2(1,1)=N2(1,1)-1; N2(end,end) = N2(end,end)-1; % subtract elements we added
              if strcmp(normaliseHist,'traj')
                N=N./N2;
              end
              if gaussianSigma<=0
                imagesc(c{1}([1 end]),c{2}([1 end]),N');     % c - pixel centres, N - pixel values (NOTE: TRANSPOSE)
              else
                imagesc(c{1}([1 end]),c{2}([1 end]),imgaussfilt(N',gaussianSigma));
              end             
              axis xy % ensure y axis points up
              colorbar
              title('Left to Right')
              xlim([ min(xyTrajTrue(:,1)), max(xyTrajTrue(:,1)) ]); ylim([ min(xyTrajTrue(:,2)), max(xyTrajTrue(:,2)) ]);
subplot(2,2,3)
[N,c] = hist3([[xyRtoL(validRLspikes,1);0;max(xyTrajTrue(:,1))], ...
              [xyRtoL(validRLspikes,2);0;max(xyTrajTrue(:,2))]], ...
              numBins ); N(1,1)=N(1,1)-1; N(end,end) = N(end,end)-1; % subtract elements we added
[N2,c2] = hist3([[xyTrajTrue(rightToLeft,1);0;max(xyTrajTrue(:,1))], ...   % Trajectoriess Histogram - i.e. how often is the trajectory in the same bin as we used for spikes
                 [xyTrajTrue(rightToLeft,2);0;max(xyTrajTrue(:,2))]], ...
                 numBins); N2(1,1)=N2(1,1)-1; N2(end,end) = N2(end,end)-1; % subtract elements we added
              if strcmp(normaliseHist,'traj')
                N=N./N2;
              end
              if gaussianSigma<=0
                imagesc(c{1}([1 end]),c{2}([1 end]),N');     % c - pixel centres, N - pixel values (NOTE: TRANSPOSE)
              else
                imagesc(c{1}([1 end]),c{2}([1 end]),imgaussfilt(N',gaussianSigma));
              end
              axis xy % ensure y axis points up
              colorbar
              title('Right to Left')
              xlim([ min(xyTrajTrue(:,1)), max(xyTrajTrue(:,1)) ]); ylim([ min(xyTrajTrue(:,2)), max(xyTrajTrue(:,2)) ]);

subplot(2,2,2)
[N,c] = hist3([[xyDtoU(validDUspikes,1);0;max(xyTrajTrue(:,1))], ...
              [xyDtoU(validDUspikes,2);0;max(xyTrajTrue(:,2))]], ...
              numBins ); N(1,1)=N(1,1)-1; N(end,end) = N(end,end)-1; % subtract elements we added
[N2,c2] = hist3([[xyTrajTrue(downToUp,1);0;max(xyTrajTrue(:,1))], ...   % Trajectoriess Histogram - i.e. how often is the trajectory in the same bin as we used for spikes
                 [xyTrajTrue(downToUp,2);0;max(xyTrajTrue(:,2))]], ...
                 numBins); N2(1,1)=N2(1,1)-1; N2(end,end) = N2(end,end)-1; % subtract elements we added
              if strcmp(normaliseHist,'traj')
                N=N./N2;
              end
              if gaussianSigma<=0
                imagesc(c{1}([1 end]),c{2}([1 end]),N');     % c - pixel centres, N - pixel values (NOTE: TRANSPOSE)
              else
                imagesc(c{1}([1 end]),c{2}([1 end]),imgaussfilt(N',gaussianSigma));
              end
              axis xy % ensure y axis points up
              colorbar
              title('Down to Up')
              xlim([ min(xyTrajTrue(:,1)), max(xyTrajTrue(:,1)) ]); ylim([ min(xyTrajTrue(:,2)), max(xyTrajTrue(:,2)) ]);
subplot(2,2,4)
[N,c] = hist3([[xyUtoD(validUDspikes,1);0;max(xyTrajTrue(:,1))], ...
              [xyUtoD(validUDspikes,2);0;max(xyTrajTrue(:,2))]], ...
              numBins ); N(1,1)=N(1,1)-1; N(end,end) = N(end,end)-1; % subtract elements we added
[N2,c2] = hist3([[xyTrajTrue(upToDown,1);0;max(xyTrajTrue(:,1))], ...   % Trajectoriess Histogram - i.e. how often is the trajectory in the same bin as we used for spikes
                 [xyTrajTrue(upToDown,2);0;max(xyTrajTrue(:,2))]], ...
                 numBins); N2(1,1)=N2(1,1)-1; N2(end,end) = N2(end,end)-1; % subtract elements we added
              if strcmp(normaliseHist,'traj')
                N=N./N2;
              end
              if gaussianSigma<=0
                imagesc(c{1}([1 end]),c{2}([1 end]),N');     % c - pixel centres, N - pixel values (NOTE: TRANSPOSE)
              else
                imagesc(c{1}([1 end]),c{2}([1 end]),imgaussfilt(N',gaussianSigma));
              end
              axis xy % ensure y axis points up
              colorbar
              title('Up to Down')
              xlim([ min(xyTrajTrue(:,1)), max(xyTrajTrue(:,1)) ]); ylim([ min(xyTrajTrue(:,2)), max(xyTrajTrue(:,2)) ]);