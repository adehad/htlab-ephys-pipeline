%% get data
fileName = '181017_06.bin';
nCh = 2;

fileID = fopen(fileName);
data = fread(fileID,[nCh, Inf], 'int16','l'); % little endian open
data = int16(data);
disp('done loading data')
size(data)

%% non filtered
cutData = data(1,763500:764000);
%cutData = data(108000:115000);
figure
%plot(cutData,'b')
hold on

dataFilt = highpass(double(cutData),300,30000,'Steepness',0.999);
plot(dataFilt,'color',[0.5 0.5 0.5])

%% filtered
dataDiscard = lowpass(double(cutData),100,30000,'Steepness',0.999);
dataFilt = cutData - dataDiscard;
plot(dataFilt,'g')

dataDiscard = lowpass(double(cutData),600,30000,'Steepness',0.999);
dataFilt = cutData - dataDiscard;
plot(dataFilt,'m')

dataFilt = highpass(double(cutData),400,30000,'Steepness',0.999);
plot(dataFilt,'y')

dataFilt = highpass(double(cutData),600,30000,'Steepness',0.999);
plot(dataFilt,'r')

legend('raw','lowpass100','lowpass600','highpass400','highpass600');