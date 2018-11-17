%% Extracts spike times to _sortedmat



extractTrialUnits();

%%
filename.raw=['C:\Users\Adehad\Desktop\dragonFlyUROP\Code\KiloSort\sample data\150526\Tetrode test data\150526__MovingObjects_1.bin'];
filename.sortOutput=['C:\Users\Adehad\Desktop\dragonFlyUROP\Code\KiloSort\sample data\150526\Tetrode test data\150526_01_sorted.mat'];


[m, fpath, mfile] = readMetafile2('150526__MovingObjects_1.meta','C:\Users\Adehad\Desktop\dragonFlyUROP\Code\KiloSort\sample data\150526\Tetrode test data\');
% [m, fpath, mfile] = readMetafile();
m.metafile = mfile;
m.metapath = fpath;
m.pdch      = m.nChans; %assume pd is last ch
m.ech       = 1:m.nChans-1; % ephys channel(s) is everything except the last
m.dbytes    = 2; % byte size of data - i.e. int16 is 2 bytes
m.msec      = m.sRateHz/1000; % conversion factor from ms time to sample number


[m,s] = extractTrialUnitWaves(filename.raw, ... % Binary File
                      filename.sortOutput, ...  % _sorted.mat file
                      m, ...                    % metafile struct, m
                      []);              % filename to store output, leave as [] if you don't want to save