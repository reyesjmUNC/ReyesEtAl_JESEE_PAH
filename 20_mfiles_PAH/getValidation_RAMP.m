function [] = getValidation_RAMP(pah)
%
% GETVALIDATION_RAMP perform a cross-validation on PAH in order to determine the
% optimum nsmax, nhmax, st metric, alphas, and weightvar 
%
% inputs: pah        - This is the name of the PAH being analyzed. This file
%                      should be in single quotes and have the same name as
%                      listed in the heading of PAHData.txt, ie 'chrysene'.
%                      Default: 'benz_a_anthracene'
%
% outputs:  
%        MSE        -   The Mean Squared Error of the X-validation.
%

% default values
if nargin < 1, pah = 1; end

% load data
load('matfiles/pah_data.mat');
idxPAH = ~isnan(val(:,pah+4));
cPM = [ProjectX ProjectY Time];

% loading covariance model to obtain baseline dmax3
load(sprintf('matfiles/covmodel_%s.mat',valname{pah+4}));
dmax3 = covparam{1,1}(2)/covparam{1,1}(3);

% loading soft data
load(sprintf('matfiles/soft_%s_RAMP.mat',valname{pah+4}));
cPM = [ProjectX ProjectY Time];

cs = cPM(~idxPAH,:);
zs = mPAHs(~idxPAH);
vs = vPAHs(~idxPAH);
zh = mPAHs(idxPAH);
nsmax = 4;

ck = cPM(idxPAH,:);
ch = cPM(idxPAH,:);

% estimation parameters
order = 0; % constant mean trend
options = BMEoptions;
options(1) = 0;
nhmax = 5;
dmax = [ 1000000 1000 dmax3 ];

% perform Xval on each unique station
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
    
    idx = ck(:,1) == unicMS(i,1) & ck(:,2) == unicMS(i,2);
    ckX{i} = ck(idx,:);
    chX{i} = ch(~idx,:);
    zhX{i}  = zh(idx,:);
    zhXval{i} = zh(~idx,:);
    
    % this is all being done is log-space
    [zkX{i},vkX{i}] = krigingME2(ckX{i},chX{i},cs,zhXval{i},zs,vs, ...
        covmodel,covparam,nhmax,nsmax,dmax,order,options);
    zkEXPX{i} = exp(zkX{i}); % take the exp of results
    
end

% saving data
save(sprintf('matfiles/BMEXVal_%s_RAMP.mat',valname{pah+4}), ...
    'ckX','chX','cs','zhX','zhXval','zs','vs','covmodel','covparam','nhmax', ...
    'nsmax','dmax','order','options','zkX','vkX','zkEXPX');

end