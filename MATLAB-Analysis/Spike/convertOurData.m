%% Adds dummy channels to Binary Files & Can: extract PD, Interlace, Invert
% Last Updated: 02/11/2018
%
% INPUT: nChans           (number of channels in binary file)
%        dataChan       (what channels are your ephys data signals)
%        adcChan          (other channel of interest, e.g. photodiode leave
%                               empty (i.e. []) if don't want to save)
%        nChansDesired    (if greater than length(dataChan) will pad with
%                           dummy channels of 1s)
%        interlaceChans   (if you want to interlace channel data)
%        invertBit        (-1 if you want to invert)
%
% OUTPUT:                 Stored is the same folder as dataFolderName
%      originalFileName_converted.bin (binary file with paddded channels
%                                    each row is a different channel)
%      originalFileName_ADC.bin        (ADC channels, never padded)


% File Name
clear fileName newFileName
    namePart1 = '150526_merged';
    namePart2 = '.bin';

% Set paramters of the Input File
    nChans = 5; % original number of channels in binary file
    dataChan = [1:4]; % which channel has your signal, array, 1 indexed - First channel is 1
    adcChan = [5];


% Set parameters of Output File
   nChansDesired = 4; % How many signal channels needed
                      % Must be larger or equal to number of channels in dataChan
                      % If larger, dummy channels are padded with 1s

   interlaceChans = 0; % 1 = yes I want to interlace
                        % by default the length(dataChan) is dividen by 2 and
                        % the second half interlaces into the first
                        % change order of interlacing by changing order of
                        % elements in dataChan array
   invertBit = 1; % 1 = NO, -1 = INVERT

% Name New File File
fileName = [ namePart1, namePart2 ];
namePart1b = '_converted';
namePart1c = '_ADC';
newFileName = [ namePart1, namePart1b, namePart2 ];
newFileNameADC = [ namePart1, namePart1c, namePart2 ];

% Open file
    clear temp
    fileID = fopen(fileName);
    temp = fread(fileID,[nChans, Inf], '*int16','l'); % little endian open
    fclose(fileID);

% Create new bin file
    clear dataRAW
    initVar = nChansDesired;
    dataRAW = ones(initVar, length(temp),'int16');

    dataRAW(1:length(dataChan),:) = temp(dataChan,:);
    dataRAW = int16(dataRAW);


% save new bin file
    fileID = fopen(newFileName,'w');
    fwrite(fileID, dataRAW, 'int16','l'); % little endian write
    fclose(fileID);

    if fileID > 0
        disp('Successfully saved SIGNAL new data');
    else
        disp('The file could not be saved, check if the newFileName is legal')
    end

    clear adcRAW
    adcRAW = zeros(length(adcChan), length(temp)); % ,'int16' % DO NOT int16
    adcRAW(1:length(adcChan),:) = temp(adcChan,:);
    adcRAW = adcRAW;    % int16() % DO NOT int16

% save new ADC bin file
    fileID = fopen(newFileNameADC,'w');
    fwrite(fileID, adcRAW, 'int16','l'); % little endian write
    fclose(fileID);

    if fileID > 0
        disp('Successfully saved ADC new data');
    else
        disp('The file could not be saved, check if the newFileNameADC is legal')
    end

clear temp