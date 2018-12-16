%% get data
fileName = '07_08_09_filtered_data_merged.bin';
fileID = fopen(fileName);
data = fread(fileID,[4, Inf], 'int16','l'); % little endian open
data = int16(data);
disp('done loading data')

%% non filtered
cutData = data(1,630000:640000);
%cutData = data(108000:115000);
figure
plot(cutData,'b')
hold on

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