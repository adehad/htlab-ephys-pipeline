uiopen('*_visual.mat');
s1 = s;
stim1 = stim;
m1 = m;
uiopen('*_visual.mat');
s2 = s;
stim2 = stim;
m2 = m;
stim.repeatIndex = [stim1.repeatIndex stim2.repeatIndex(2:end)+stim1.repeatIndex(end)];
stim.stimLength = [stim1.stimLength stim2.stimLength];
stim.stimLength = [stim1.stimLength stim2.stimLength];