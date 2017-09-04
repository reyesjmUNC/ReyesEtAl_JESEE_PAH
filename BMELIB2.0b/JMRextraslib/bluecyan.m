function h = bluecyan(m)
%BLUECYAN    Black-blue-cyan-white color map, or "deep blue"
%   BLUECYAN(M) returns an M-by-3 matrix containing a "bluecyan" colormap.
%   BLUECYAN, by itself, is the same length as the current colormap.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(bluecyan)
%
%   See also HSV, GRAY, PINK, COOL, BONE, COPPER, FLAG, COLORMAP, RGBPLOT,
%   REDYELLOW, REDPINK, GREENYELLOW, GREENCYAN, BLUEPINK 


%   M. Serre, 3-01-07,

if nargin < 1, m = size(get(gcf,'colormap'),1); end
rgb=flipud(hot(m));
h = rgb(:,[3 2 1]);
