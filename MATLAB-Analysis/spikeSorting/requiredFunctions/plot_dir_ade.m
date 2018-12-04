function [h1, h2] = plot_dir_ade (vX, vY,varargin)
%function [h1, h2] = plot_dir (vX, vY)
%Plotting x-y variables with direction indicating vector to the next element.
%Example
%   vX = linspace(0,2*pi, 10)';
%   vY = sin (vX);
%   plot_dir(vX, vY);
% varargin{1} = arrow head scaling default 0.5
% varargin{2} = axis handle

% disp(length(varargin));
if length(varargin)==2
    if ~isempty(varargin{1})
        rMag = varargin{1};
    else
        rMag = 0.5;
    end
    if ~isempty(varargin{2})
        ax = varargin{2};
    end
elseif length(varargin)==1
    if ~isempty(varargin{1})
        rMag = varargin{1};
    else
        rMag = 0.5;
    end
elseif ~isempty(varargin)
    rMag = 0.5;
else
	rMag = 0.5;
end

% Length of vector
lenTime = length(vX);

% Indices of tails of arrows
vSelect0 = 1:(lenTime-1);
% Indices of tails of arrows
vSelect1 = vSelect0 + 1;

% X coordinates of tails of arrows
vXQ0 = vX(vSelect0, 1);
% Y coordinates of tails of arrows
vYQ0 = vY(vSelect0, 1);

% X coordinates of heads of arrows
vXQ1 = vX(vSelect1, 1);
% Y coordinates of heads of arrows
vYQ1 = vY(vSelect1, 1);

% vector difference between heads & tails
vPx = (vXQ1 - vXQ0) * rMag;
vPy = (vYQ1 - vYQ0) * rMag;


if exist('ax')
   % make plot 
%     h1 = plot (ax,vX, vY, '.-');
    % add arrows 
    h2 = quiver (ax,vXQ0,vYQ0, vPx, vPy, 0, 'b'); 
else
    % make plot 
    h1 = plot (vX, vY, '.-');
    % add arrows 
    h2 = quiver (vXQ0,vYQ0, vPx, vPy, 0, 'r');
end
