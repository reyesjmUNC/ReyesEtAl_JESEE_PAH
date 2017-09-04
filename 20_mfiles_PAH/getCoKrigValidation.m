function [] = getCoKrigValidation(pah)
% GETVALIDATION This function will perfrom a cross validation on PAH using
% cokriging
%
% inputs: pah        - This is the name of the PAH being analyzed. This file
%                      should be in single quotes and have the same name as
%                      listed in the heading of PAHData.txt, ie 'chrysene'.
%                      Default: 'benz_a_anthracene'
%                      Type: string

% default values
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

% perform Xval on each unique station
ck = chPM(idx,:);
ch = chPM(idx,:);
unicMS = unique(ck(:,1:2),'rows');
len = size(unicMS,1);
ckX = cell(len,1);
chX = cell(len,1);
zhX = cell(len,1);
zhXval = cell(len,1);
zkX = cell(len,1);
vkX = cell(len,1);
zkEXPX = cell(len,1);
for i = 1:len
    
    idx = chjustPAH(:,1) == unicMS(i,1) & chjustPAH(:,2) == unicMS(i,2);
    zhX{i}  = zjustPAH(idx,:);    
    chPAHtemp = chjustPAH(idx,:);
    ckPAHtemp = chjustPAH(~idx,:);
    zPAHtemp = zjustPAH(~idx,:);
    
    ckX{i} = { [chPAHtemp] [1*ones(size(chPAHtemp,1),1)] }; 
    chX{i} = { [ckPAHtemp;chPM] [1*ones(size(ckPAHtemp,1),1);2*ones(length(chPM),1)] };
    cs = [];    
    zhXval{i} = [ zPAHtemp ; zPM ];
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
    [moments,info]=BMEprobaMoments(ckX{i},chX{i},cs,zhXval{i}, ...
        softpdftype,nl,limi,probdens,covmodel,covparam,nhmax,nsmax,dmax,order,options);
    zkX{i} = moments(:,1);
    vkX{i} = moments(:,2);
    zkEXPX{i} = exp(zkX{i}); % take the exp of results

end

% saving data
save(sprintf('matfiles/BMEXVal_%s_cokrig.mat',valname{pah+4}), ...
    'ckX','chX','cs','zhX','zhXval','nl','limi','probdens','covmodel', ...
    'covparam','nhmax','nsmax','dmax','order','options','zkX','vkX','zkEXPX');

end