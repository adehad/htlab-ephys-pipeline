%% Converting Open-Ephs Data (organised in folders) to Binary Files needed for Spike Sorting

% Folder separator: Windows = '\' Linux = '/'
    pathToDataFolder    = '\\bg-thefarm-2012\Shares\Lin_group_data\DATA\Ephys\testStimulus\';   % If script already can reach dataFolder, set to ''
    dataFolderNames     = 'all';% String array (DOUBLE QUOTES) such as ["181108",...] OR 'all' 
    overwriteFiles      = 1;    % 1 = overwrite, 0 = do not overwrite
    
% Data Channel(s) - 100_CH<>.continuous files
    dataCh              = [];  % Array
% PhotoDiode Channel - 100_ADC<>.continuous files
    adcCh               = 1;    % Array
    
% Output parameters
    nChDesired          = 4;    % Number of signal channels needed
                                % Must be larger or equal to number of channels in dataChan
                                % If larger, dummy channels are padded with 1s
    interlaceCh         = 0;    % 1 = interlace, 0 = do not interlace
                                % by default the length(dataChan) is divided by 2 and
                                % the second half interlaces into the first
                                % change order of interlacing by changing order of
                                % elements in dataChan array
    invertCh            = 0;    % 1 = invert, 0 = do not invert
    mergeCh             = 1;    % 1 = merge data, 0 = do not merge data

OEtoBin(pathToDataFolder,dataFolderNames,overwriteFiles,dataCh,adcCh,nChDesired,interlaceCh,invertCh,mergeCh)

%% Converting .bin files (organised in a single folder) to Binary Files needed for Spike Sorting

% Folder separator: Windows = '\' Linux = '/'
    pathToDataFolder    = 'E:\Box Sync\DragonVision_DanielKo\Data\Ephys\testEnv';   % If script already can reach dataFolder, set to ''
    dataNames     = 'all';% String array (DOUBLE QUOTES) such as ["181108",...] OR 'all' 
    overwriteFiles      = 1;    % 1 = overwrite, 0 = do not overwrite
    
% Data Channel(s) - 100_CH<>.continuous files
    dataCh              = 4;  % Number of channels in input .bin
    
% Output parameters
    nChDesired          = 4;    % Number of signal channels needed
                                % Must be larger or equal to number of channels in dataChan
                                % If larger, dummy channels are padded with 1s
    interlaceCh         = 0;    % 1 = interlace, 0 = do not interlace
                                % by default the length(dataChan) is divided by 2 and
                                % the second half interlaces into the first
                                % change order of interlacing by changing order of
                                % elements in dataChan array
    invertCh            = 0;    % 1 = invert, 0 = do not invert
    mergeCh             = 1;    % 1 = merge data, 0 = do not merge data

binToBin(pathToDataFolder,dataNames,overwriteFiles,dataCh,nChDesired,interlaceCh,invertCh,mergeCh)