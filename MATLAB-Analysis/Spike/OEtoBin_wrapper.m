%% Converting Open-Ephs Data (organised in folders) to Binary Files needed for Spike Sorting
clear pathToDataFolder dataFolderName
% Change Folder separator from '\' to '/' if running on Linux
% You MUST end the folder names with the folder separator    
    pathToDataFolder = '';  % If script already can reach dataFolder, set to ''
    dataFolderName = ['181108',]; 
        % Array
        
% Set paramters of data
% Data Channel(s)   - 100_CH<>.continuous files
    dataChan = [];     % Array
% PhotoDiode Channel    - 100_ADC<>.continuous files
    adcChan = [1];     % Array

% Set parameters of output
    nChansDesired = 1; % How many signal channels needed
                       % Must be larger or equal to number of channels in dataChan
                       % If larger, dummy channels are padded with 1s
    interlaceChans = 0; % 1 = yes I want to interlace
                        % by default the length(dataChan) is dividen by 2 and
                        % the second half interlaces into the first
                        % change order of interlacing by changing order of
                        % elements in dataChan array
    invertChans = 1; % 1 = NO, -1 = INVERT

    concatCh = 1; % 1 = Yes concatenate, anything else = No

OEtoBin_wrapper();