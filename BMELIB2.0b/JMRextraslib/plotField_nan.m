function [xg yg Zg] = plotField(sk,zk,ax,maskcontour,colorz)
% plotField       - Makes a color map of the field of values
%
% SYNTAX:
%   plotField(sk,zk,ax,maskcontour)
%
% INPUT
%   sk            nk x 2      matrix of estimation points
%   zk            nk x 1      column vector of field values at the points sk
%   ax            1 x 4       [xmin xmax ymin ymax] of the area over which to plot the field
%                             default value is ax=[], which will use an area over all the field
%   maskcontour   n x 2       matrix of points defining the outer contour of the mask
%                             default value is maskcontour=[], which will not use any mask


if nargin < 3, ax=[]; end;            
if nargin < 4, maskcontour=[]; end;
if nargin < 5, colorz = colormap(jet); end

nxpix=148;           % Number of pixels in the x-direction used to create the color map 
nypix=112;           % Number of pixels in the y-direction used to create the color map 
maskfillcolor='w';   % character defining the color to use to fill outside of the mask
                     % 'w' is for white, see help plot for other colors 
masklinetype='k';    % character defining the color to use for the linetype of 
                     % the mask contour. 'k' is for black
                     

if isempty(ax)
  ax=[min(sk(:,1))-2 max(sk(:,1))+2 min(sk(:,2))-2 max(sk(:,2))+2];
end;

dx=diff(ax(1:2))/nxpix/2;
dy=diff(ax(3:4))/nypix/2;
ax=[ax(1)+dx ax(2)-dx ax(3)+dy ax(4)-dy]; 


dx1=diff(ax(1:2))/nxpix;
dy1=diff(ax(3:4))/nypix;
% xg=[ax(1):dx1:ax(2)+dx1-eps];
% yg=[ax(3):dy1:ax(4)+dy1-eps];
% [xg yg]=meshgrid(xg,yg);                   % Gridpoint of pixel used to display the field
% Zg=griddata(sk(:,1),sk(:,2),zk,xg,yg);     % Value of the field at the pixel gridppoints
% Zg=reshape(Zg,size(xg));
% maxZg=max(max(Zg));

xg = unique(round(sk(:,1)));
yg = unique(round(sk(:,2)));
[xg yg] = meshgrid(xg,yg);
Zg = reshape(zk,size(xg));
maxZg = max(max(Zg));

figure;
pcolor(xg,yg,Zg);        % Create the color map
colormap(colorz);
%shading interp;
shading flat;
hold on

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
