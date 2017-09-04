function [] = getValidation(pah,LR,soft)
%
% GETVALIDATION perform a cross-validation on PAH in order to determine the
% optimum nsmax, nhmax, st metric, alphas, and weightvar 
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
%        MSE        -   The Mean Squared Error of the X-validation.
%

% default values
if nargin < 1, pah = 1; end
if nargin < 2, LR = 1; end 
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
save(sprintf('matfiles/BMEXVal_%s_ellip_%s_%s.mat',valname{pah+4},LRstr,hardstr), ...
    'ckX','chX','cs','zhX','zhXval','zs','vs','covmodel','covparam','nhmax', ...
    'nsmax','dmax','order','options','zkX','vkX','zkEXPX');

end