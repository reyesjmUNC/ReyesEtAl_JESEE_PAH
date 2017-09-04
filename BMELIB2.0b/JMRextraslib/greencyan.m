function h = greencyan(m)
%GREENCYAN    Black-green-cyan-white color map, or "algea"
%   GREENCYAN(M) returns an M-by-3 matrix containing a "greencyan" colormap.
%   GREENCYAN, by itself, is the same length as the current colormap.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(greencyan)
%
%   See also HOT, HSV, GRAY, PINK, COOL, BONE, COPPER, FLAG, COLORMAP, RGBPLOT,
%   REDYELLOW, REDPINK, GREENYELLOW, BLUECYAN, BLUEPINK 


%   M. Serre, 3-01-07,

if nargin < 1, m = size(get(gcf,'colormap'),1); end
rgb=flipud(hot(m));
h = rgb(:,[3 1 2]);
