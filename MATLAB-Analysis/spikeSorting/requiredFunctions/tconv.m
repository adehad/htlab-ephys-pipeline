function [dnew,fmax] = tconv(data,kernel,cmid)
%
% conv data with kernel and compensate for time-shift
% induced by convolution
%
% length(dnew) = length(data)
%
% assume kernel center is at max unless cmid=1
%
%IN: data = time series to be convolved
%   kernel = convolution kernel
%   cmid = assume middle of kernel is center (ie for derivative)
%
%OUT: dnew = convolved, time-shifted data
%     fmax = center of kernel
%
% AL, Caltech, 8/00
%
% [dnew,fmax] = tconv(data,kernel,cmid)
%

if ~exist('cmid','var')  || isempty(cmid),  cmid = 0; end;
dnew = conv(data,kernel);
if (cmid == 0)
    fmax = find(kernel == max(kernel),1);
else fmax = round(length(kernel)/2);
end;
dnew = dnew(fmax:fmax+length(data)-1);

% recenter even-length (non-symmetric) filter with cmid centering)
if (mod(length(kernel),2) == 0) & cmid == 1
    dnewi = interp(dnew,2);
    dnew = dnewi(2:2:end);
end;
