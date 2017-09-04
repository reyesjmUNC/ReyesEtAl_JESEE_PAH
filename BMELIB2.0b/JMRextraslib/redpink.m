function h = redpink(m)
%REDPINK    Black-red-pink-white color map, or "blood"
%   REDPINK(M) returns an M-by-3 matrix containing a "redpink" colormap.
%   REDPINK, by itself, is the same length as the current colormap.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(redpink)
%
%   See also HOT, HSV, GRAY, PINK, COOL, BONE, COPPER, FLAG, COLORMAP, RGBPLOT,
%   REDYELLOW, GREENYELLOW, GREENCYAN, BLUECYAN, BLUEPINK 


%   M. Serre, 3-01-07,

if nargin < 1, m = size(get(gcf,'colormap'),1); end
rgb=flipud(hot(m));
h = rgb(:,[1 3 2]);
