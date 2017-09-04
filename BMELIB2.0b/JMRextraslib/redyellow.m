function h = redyellow(m)
%REDYELLOW    Black-red-yellow-white color map, or "hot reversed"
%   REDYELLOW(M) returns an M-by-3 matrix containing a "redyellow" colormap.
%   REDYELLOW, by itself, is the same length as the current colormap.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(redyellow)
%
%   See also HOT, HSV, GRAY, PINK, COOL, BONE, COPPER, FLAG, COLORMAP, RGBPLOT,
%   REDPINK, GREENYELLOW, GREENCYAN, BLUECYAN, BLUEPINK 


%   M. Serre, 3-01-07,

if nargin < 1, m = size(get(gcf,'colormap'),1); end
h=flipud(hot(m));
