function [metafile, PathName,FileName] = readMetafile

%select desired metafile from folder
[FileName,PathName] = uigetfile('*.meta','Select any metafile');
y = [PathName,FileName];

%open selected file, read data using string format, delimeter is a space
fid = fopen(strcat(PathName,FileName));
a = textscan(fid,'%s','delimiter','','whitespace',''); %,'BufSize',8190);

%open cell "a"(1x1), "a" contains an array of the parameters from
%metafile
array = a{1,1};


%loop through the entire array 
%go through each row, convert value in that cell into a string, find the
%index of the "=" in each row, field is text before "=" while value is
%string after "="

metafile=struct;
for i=1:length(array)
    str = array(i,1);
    str = str{1,1};
    
    equalind = regexpi(str, '=');
    fieldvalue = str(equalind+2:end);
    fieldvalue_num = str2double(fieldvalue);
    if ~isnan(fieldvalue_num)
        fieldvalue = fieldvalue_num;
    end;
    field = str(1:equalind-1);
    metafile = setfield(metafile,field,fieldvalue);  
    
end

fclose(fid);

end



%Note:must use str2num to get numerical values back into original form,
%when calling certain field values
%metafile=ans;
%str2num(Metafile.fileTimeSecs)
