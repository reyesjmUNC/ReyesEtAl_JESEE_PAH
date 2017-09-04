function [] = getBMEest_cokrig(pah)
% GETBMEEST_COKRIG performs the cokriging estimates on PAH data. Note that for PAH, soft 
% data is only defined where PM2.5 is measured. However, for graphing
% purposes, estimation of PAH is needed on a regular grid.
%
% inputs: pah        - This is the name of the PAH being analyzed. This file
%                      should be in single quotes and have the same name as
%                      listed in the heading of PAHData.txt, ie 'chrysene'.
%                      Default: 'benz_a_anthracene'
%                      Type: string
%
% outputs:  
% There are no outputs.
%

if nargin < 1, pah = 1; end

% load data
load('matfiles/pah_data.mat');
zPAH = val(:,pah+4);
idx = ~isnan(zPAH);
zPAH(~idx) = NaN; zPAH(idx) = log(zPAH(idx)); 
zjustPAH = zPAH; zjustPAH(~idx) = [];
zPM = log(val(:,4));
chPAH = [ProjectX ProjectY Time]; chPAH(~idx,:) = NaN; 
chjustPAH = chPAH; chjustPAH(~idx,:) = [];
chPM = [ProjectX ProjectY Time];

% loading covariance model
load(sprintf('matfiles/covmodel_cokrig%s.mat',valname{pah+4}));
dmax3 = covparam{1}{2}(1)/covparam{1}{2}(2);

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
idx = ~isnan(val(:,pah+4));
DaysSince = unique(Time(idx));

len = length(DaysSince);
temp = repmat(DaysSince,length(xg(:)),1);
ck=[repmat(xg(:),len,1) repmat(yg(:),len,1) temp(:)];

% all cokriging parameters
ck = { [ck] [1*ones(size(ck,1),1)] }; 
ch = { [chjustPAH;chPM] [1*ones(size(chjustPAH,1),1);2*ones(length(chPM),1)] };  
cs = [];  
zh = [ zjustPAH ; zPM ];
softpdftype = 1;
nl = [];
limi = [];
probdens = [];
nhmax = [ 5 ; 5 ];
nsmax = [0;0];
dmax = [ 1000000 600 dmax3 ];
order=[0;0];
options = BMEoptions;
options(1) = 0;
options(3) = 150000;

% this is all being done is log-space
% tic/toc ~= 140 seconds
[moments,info]=BMEprobaMoments(ck,ch,cs,zh,softpdftype,nl,limi, ...
    probdens,covmodel,covparam,nhmax,nsmax,dmax,order,options);
zk = moments(:,1);
vk = moments(:,2);
zkEXP = exp(zk); % take the exp of results

% saving results
save(sprintf('matfiles/BMEest_%s_cokrig.mat',valname{pah+4}), ...
    'ck','ch','cs','zh','softpdftype','nl','limi', 'probdens', ...
    'covmodel','covparam','nhmax','nsmax','dmax','order','options','zk','vk','zkEXP');

end