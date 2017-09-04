function h = greenyellow(m)
%GREENYELLOW    Black-green-yellow-white color map.
%   GREENYELLOW(M) returns an M-by-3 matrix containing a "greenyellow" colormap.
%   GREENYELLOW, by itself, is the same length as the current colormap.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(greenyellow)
%
%   See also HOT, HSV, GRAY, PINK, COOL, BONE, COPPER, FLAG, COLORMAP, RGBPLOT,
%   REDYELLOW, REDPINK, GREENCYAN, BLUECYAN, BLUEPINK 


%   M. Serre, 3-01-07,

if nargin < 1, m = size(get(gcf,'colormap'),1); end
rgb=flipud(hot(m));
h = rgb(:,[2 1 3]);
