function [] = getBMEest(pah,LR,soft)
% GETBMEEST performs the BME estimates on PAH data. Note that for PAH, soft 
% data is only defined where PM2.5 is measured. However, for graphing
% purposes, estimation of PAH is needed on a regular grid.
%
% inputs: pah        - This is the name of the PAH being analyzed. This file
%                      should be in single quotes and have the same name as
%                      listed in the heading of PAHData.txt, ie 'chrysene'.
%                      Default: 'benz_a_anthracene'
%                      Type: string
%         LR         - This variable determines which method is used to
%                      find the soft data. If LR = 1, a linear regression 
%                      method will be used to determine the soft data. If 
%                      LR = 0, a mass fraction apprach will be used to 
%                      determine the soft data.
%         soft       - This variable determines hard or soft data used. 
%                      soft = 1 means soft and soft = 0 means hard.
%
% outputs:  
% There are no outputs.
%

if nargin < 1, pah = 1; end
if nargin < 2, LR = 0; end
if nargin < 3, soft = 1; end

% load data
load('matfiles/pah_data.mat');
idxPAH = ~isnan(val(:,pah+4));
cPM = [ProjectX ProjectY Time];

% loading covariance model to obtain baseline dmax3
load(sprintf('matfiles/covmodel_%s.mat',valname{pah+4}));
dmax3 = covparam{1,1}(2)/covparam{1,1}(3);

if soft == 0 % no soft data
    
    cs = [];
    zs = [];
    vs = [];
    zh = log(val(idxPAH,pah+4)); 
    nsmax = 0;
    hardstr = 'hard';
    LRstr = 'NA';
    
    
elseif soft == 1 % soft data
    
    % set strings
    if LR == 1, LRstr = 'LR'; else LRstr = 'MF'; end
    hardstr = 'soft';
    
    % loading soft data
    load(sprintf('matfiles/soft_%s_ellip_%s.mat',valname{pah+4},LRstr));
    cPM = [ProjectX ProjectY Time];
    
    cs = cPM(~idxPAH,:);
    zs = mPAHs(~idxPAH);
    vs = vPAHs(~idxPAH);
    zh = mPAHs(idxPAH);
    nsmax = 4;
  
end

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

ch = cPM(idxPAH,:);

% estimation parameters
order = 0; % constant mean trend
options = BMEoptions;
options(1) = 0;
nhmax = 5;
dmax = [ 1000000 1000 dmax3 ];

% performing the BME estimates
% tic/toc ~= 130 seconds
[zk,vk] = krigingME2(ck,ch,cs,zh,zs,vs,covmodel,covparam,nhmax,nsmax,dmax,order,options);
zkEXP = exp(zk);

% saving results
save(sprintf('matfiles/BMEest_%s_ellip_%s_%s.mat',valname{pah+4},LRstr,hardstr), ...
    'ck','ch','cs','zh','zs','vs', 'covmodel','covparam','nhmax', ...
    'nsmax','dmax','order','options','zk','vk','zkEXP');

end