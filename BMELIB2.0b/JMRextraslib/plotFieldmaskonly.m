function [] = plotFieldmaskonly(ax,maskcontour)
% plotField       - Makes a mask only
%
% SYNTAX:
%   plotField(ax,maskcontour)
%
% INPUT
%   ax            1 x 4       [xmin xmax ymin ymax] of the area over which to plot the field
%                             default value is ax=[], which will use an area over all the field
%   maskcontour   n x 2       matrix of points defining the outer contour of the mask
%                             default value is maskcontour=[], which will not use any mask


if nargin<1, ax=[]; end;            
if nargin<2, maskcontour=[]; end;

maskfillcolor='w';   % character defining the color to use to fill outside of the mask
                     % 'w' is for white, see help plot for other colors 
masklinetype='k';    % character defining the color to use for the linetype of 
                     % the mask contour. 'k' is for black
                     
%  Create a mask and fill out the outside of the mask with a uniform color
if ~isempty(maskcontour)
    axmask = [ min([maskcontour(:,1);ax(1)]) max([maskcontour(:,1);ax(2)])...
        min([maskcontour(:,2);ax(3)]) max([maskcontour(:,2);ax(4)]) ];
  
    x = maskcontour(:,1);
    y = maskcontour(:,2);
    xi=axmask(1); xe=axmask(2);
    yi=axmask(3); ye=axmask(4);

    i=find(x==min(x)); i=i(1);

    x=x(:);
    y=y(:);
    x=[x(i:end)' x(1:i-1)' x(i)];
    y=[y(i:end)' y(1:i-1)' y(i)];

    x=[xi   xi xe xe xi xi   x(1) x];
    y=[y(1) ye ye yi yi y(1) y(1) y];

    h=fill(x,y,maskfillcolor);
    set(h,'edgecolor','none');
    plot(maskcontour(:,1),maskcontour(:,2),masklinetype);
end
