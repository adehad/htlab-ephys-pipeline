function chans = getKwikData(kwikFile)
igroup=0; basename = 'XXX';
spkTimes = hdf5read(kwikFile, ['/channel_groups/' num2str(igroup) '/spikes/time_samples']);
spkClus = hdf5read(kwikFile, ['/channel_groups/' num2str(igroup) '/spikes/clusters/main']);
cellIDs = unique(spkClus);
ncells  = length(cellIDs);

% spkClus(spkTimes>endTime) = [];
% spkTimes(spkTimes>endTime) = [];
% spkTimes = double(spkTimes);
% spkTimes = spkTimes - startTime;
% spkClus(spkTimes<0) = [];
% spkTimes(spkTimes<0) = [];

temp =  h5info(kwikFile,  '/recordings/0/');
sampleRate = double((temp.Attributes(3).Value));
spkTimes = double(spkTimes)./sampleRate;
noise_list = [];
for icell = 1:ncells
    chans(icell).spiketimes = spkTimes(spkClus==cellIDs(icell));
%     chans(icell).ichan = igroup;
%     chans(icell).iexp = iexp;
    chans(icell).icell = cellIDs(icell);
    chans(icell).sampleRate = sampleRate;
    
    temp = h5info(kwikFile, ['/channel_groups/' num2str(igroup) '/clusters/main/' num2str(cellIDs(icell)) '/']);
    
    for idx = 1:length(temp.Attributes)
        if strcmp(temp.Attributes(idx).Name, 'cluster_group')
            ilabel = temp.Attributes(idx).Value;
            break
        end
        idx = idx + 1;
    end
    if idx <= length(temp.Attributes)
        ilabel = temp.Attributes(idx).Value;
        temp = h5info(kwikFile, ['/channel_groups/' num2str(igroup) '/cluster_groups/main']);
        labelType = eval(['temp.Groups(' num2str(ilabel+1) ').Attributes(end).Value']);
        if strcmp(labelType,'Noise')
            noise_list = [noise_list icell];
        end
        chans(icell).id = [labelType '_cluster' num2str(cellIDs(icell))];
    else
        chans(icell).id = ['none_cluster' num2str(cellIDs(icell))];
    end
end
chans(noise_list) = [];
