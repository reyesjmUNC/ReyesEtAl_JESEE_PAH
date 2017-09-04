function h = bluepink(m)
%BLUEPINK    Black-blue-pink-white color map.
%   BLUEPINK(M) returns an M-by-3 matrix containing a "bluepink" colormap.
%   BLUEPINK, by itself, is the same length as the current colormap.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(bluepink)
%
%   See also HOT, HSV, GRAY, PINK, COOL, BONE, COPPER, FLAG, COLORMAP, RGBPLOT,
%   REDYELLOW, REDPINK, GREENYELLOW, GREENCYAN, BLUECYAN 


%   M. Serre, 3-01-07,

if nargin < 1, m = size(get(gcf,'colormap'),1); end
rgb=flipud(hot(m));
h = rgb(:,[2 3 1]);
