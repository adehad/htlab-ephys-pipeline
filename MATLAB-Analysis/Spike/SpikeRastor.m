%%%% SpikeRastor
%%%% Huai-Ti Lin [Dec 2013]
%%%% This script produces a rastor plot from sorted data
%%%% INPUT: m.waves m.spikes m.xytraj m.pd m.StimGL_nloops

%% REMEMBER to import correct the StartFrame variable from the stimulus mat file


%% Basic info
D=34; %55; Distance of the DF head to the screen
W=125; %123; Width of the projection
H=70;  % Height of the projection
xMap=1280/W; %111; % pix/mm 125mm
yMap=720/H;  %63; % pix/mm % 80mm
dt=1/360;
stimLength=size(m.xytraj,1);
sort_bit=0; %to plot raster and RF for sorted and non-sorted data
% waves=m.waves;
% spikes=m.spikes;
waves=double(m.waves);
% spikes=double(m.spikes);
spikes = double(m.spikes - m.pd(1) + 1); %Start aligned the spikes  /m.msec;
pd_diff = double(diff(m.pd));

Repeat_thrs = 0.5; % this is the cutoff threshold for how consistent the spikes are over repeated trials. 

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
[RepeatSpace, RepeatIndex]=find(pd_diff>2000); % this is in 360Hz unit because diff is actually in 360Hz 
nloops = m.StimGL_nloops;
xytraj =m.xytraj;
% StartFs =StartFrame;
TarTraj_s=TarTraj;
TarTraj_dot_s=TarTraj_dot;
for nl=1:nloops-1;
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
if sort_bit==0;
spikes_360=round(spikes*360/40000);
elseif sort_bit==1;
spikes_360=round(spikes(UnitInd)*360/40000);  
end
RasterMask=zeros(1,stimLength);
ind=find(spikes_360>0 & spikes_360<stimLength);
RasterMask(spikes_360(ind))=ones(1,size(ind,2));

Raster(1,:)=RasterMask; 
SpikeTiming=find(Raster(1,:));
T=RasterMask; % linear raster array
for re=1:nloops-1;
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
subplot(10,1,1:2); plot(TarTraj(:,1),'r'); hold on; plot(TarTraj(:,2),'b');
xlim([0 stimLength])
ylim([-60 60])
subplot(10,1,3:4); plot(TarTraj_dot(:,1),'r'); hold on; plot(TarTraj_dot(:,2),'b');
xlim([0 stimLength])
ylim([-150 150])

subplot(10,1,5:9); rasterplot(TE,nloops-1,stimLength,gca);
% rasterplot(TE,nloops-1,stimLength,'plotwidth',0.5);
subplot(10,1,10); 
%bar(centers./360,nelements);
plot((1:stimLength)./360,SpikeHis)
xlim([0 stimLength./360])
ylim([0 50])
