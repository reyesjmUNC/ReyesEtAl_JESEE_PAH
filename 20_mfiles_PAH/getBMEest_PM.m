function [] = getBMEest_PM()
% GETBMEEST_PM performs the BME estimates on PM2.5 data. For graphing
% purposes, estimation of PM is needed on a regular grid.

% load data/covariance model
load('matfiles/pah_data.mat');
load('matfiles/covmodel_PM2p5.mat');

% determining the space/time locations of ck
ax = [ -85 -75 33.5 36.6 ]; % these are hard coded, specific for NC
axx = [-85;-75]; axy = [33.5;36.6];
% project longtiude/latitude
cd ../09_mfiles_projections
projxy = ell2lambertcc([axx,axy],'whiproj2001');
cd ../20_mfiles_PAH
ax = [projxy(1,1) projxy(2,1) projxy(1,2) projxy(2,2)];
nx = 25;
ny = 10;
[xg yg]=meshgrid(ax(1):diff(ax(1:2))/nx:ax(2),ax(3):diff(ax(3:4))/ny:ax(4));

% getting the estimation times, every day where PAH data exist
idx = ~isnan(val(:,5));
DaysSince = unique(Time(idx));

len = length(DaysSince);
temp = repmat(DaysSince,length(xg(:)),1);
ck=[repmat(xg(:),len,1) repmat(yg(:),len,1) temp(:)];

dmax3 = ( (f.alp*f.ar1)+((1-f.alp)*f.ar2) ) / ( (f.alp*f.at1)+((1-f.alp)*f.at2) );
dmax = [1000000 1000 dmax3];

ch = [ProjectX ProjectY Time];
cs = [];
zh = log(val(:,4));
softpdftype = 1; 
zs = [];
vs = [];
nhmax = 5;
nsmax = 0;
order = 0;
options = BMEoptions;
options(1) = 0;
options(3) = 150000;

% performing the BME estimates
% tic/toc ~= 110 seconds
[zk,vk] = krigingME2(ck,ch,cs,zh,zs,vs,covmodel,covparam,nhmax,nsmax,dmax,order,options);
zkEXP = exp(zk);

% saving results
save('matfiles/BMEest_krig_PM2p5.mat','ck','ch','cs','zh','zs','vs', ...
    'covmodel','covparam','nhmax','nsmax','dmax','order','options','zk','vk','zkEXP');

end