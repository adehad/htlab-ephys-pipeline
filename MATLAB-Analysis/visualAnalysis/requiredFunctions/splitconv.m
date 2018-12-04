function y = splitconv(d,f)
%
% convolution with minimized edge artifacts
% assumes length(d) >> length(f)
%
%IN: d = time series (vector)
%    f = kernel
%
%OUT: y = filtered d
%
% AL, janelia; 11/09
%
%y = splitconv(d,f)
%


dmid = round(size(d,1)/2);
da = d - d(1);
db = d - d(end);
daf = tconv(da,f) + d(1);
dbf = tconv(db,f) + d(end);
daf = daf(:);
dbf = dbf(:);
y = [daf(1:dmid)' dbf(dmid+1:end)'];


