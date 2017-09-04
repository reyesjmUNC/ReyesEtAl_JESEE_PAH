function [] = exhaustiveSoftDevelopment(pahnum)
% Runs 'callRegression.m' for different values of spatrad, temprad, nmax, 
% and stmetric. This is the analysis function for the LR and MF models. 
% This function will capture the MSE for LR, MSE for MF, AIC for LR, and 
% AIC for MF.

if nargin < 1, pahnum = 1; end

% load data
load('matfiles/pah_data.mat');
Time  = Time;
ProjectX = ProjectX;
ProjectY = ProjectY;

% matlabpool open 12 % parallel computing

% initialize variables
spatrad = 10000:10000:1000000; % units: meters
temprad = 3:3:350; % units: days
nmax = 1:84;
stmetric = 0.01:0.01:.9;
stmetric = 1:1:900; % changing s/t metric when going from degrees to meters
stmetric = 1000:1000:900000; % 1km:1km:900km

% loop through cylindrical/ellipsoidal neighborhood
for j = 1:2

    if j == 1 % for cylindrical neighborhood
            MSELR = NaN*ones(length(spatrad),length(temprad));
            lbMSELR = NaN*ones(length(spatrad),length(temprad));
            ubMSELR = NaN*ones(length(spatrad),length(temprad));
            AICLR = NaN*ones(length(spatrad),length(temprad));
            MSEMF = NaN*ones(length(spatrad),length(temprad));
            lbMSEMF = NaN*ones(length(spatrad),length(temprad));
            ubMSEMF = NaN*ones(length(spatrad),length(temprad));
            AICMF = NaN*ones(length(spatrad),length(temprad));

            for k = 1:length(spatrad)
                disp([j k]);
                tic
                parfor l = 1:length(temprad) 
                    [MSELR(k,l),AICLR(k,l),MSEMF(k,l),AICMF(k,l), ...
                        lbMSELR(k,l),ubMSELR(k,l),lbMSEMF(k,l),ubMSEMF(k,l)] = ...
                        calcRegression([ProjectX ProjectY Time],val(:,4),val(:,pahnum+4), ...
                        spatrad(k),temprad(l),max(nmax),max(stmetric));
                end
                toc
            end
            save(sprintf('matfiles/exhaustive_%s_cylin.mat',valname{pahnum+4}),'spatrad', ...
                'temprad','MSELR','AICLR','MSEMF','AICMF','lbMSELR','ubMSELR', ...
                'lbMSEMF','ubMSEMF');      
    else % for ellipsoidal neighborhood 
        MSELR = NaN*ones(length(nmax),length(stmetric));
        lbMSELR = NaN*ones(length(nmax),length(stmetric));
        ubMSELR = NaN*ones(length(nmax),length(stmetric));
        AICLR = NaN*ones(length(nmax),length(stmetric));
        MSEMF = NaN*ones(length(nmax),length(stmetric));
        lbMSEMF = NaN*ones(length(nmax),length(stmetric));
        ubMSEMF = NaN*ones(length(nmax),length(stmetric));
        AICMF = NaN*ones(length(nmax),length(stmetric));

        for k = 1:length(nmax)
            disp([j k]);
            tic
            for l = 1:length(stmetric) % changed from 'parfor'
                [MSELR(k,l),AICLR(k,l),MSEMF(k,l),AICMF(k,l), ...
                    lbMSELR(k,l),ubMSELR(k,l),lbMSEMF(k,l),ubMSEMF(k,l)] = ...
                    calcRegression([ProjectX ProjectY Time],val(:,4),val(:,pahnum+4), ...
                    max(spatrad),max(temprad),nmax(k),stmetric(l));
            end
            toc
        end
        save(sprintf('matfiles/exhaustive_%s_ellip_newstmetric.mat',valname{pahnum+4}),'nmax', ...
            'stmetric','MSELR','AICLR','MSEMF','AICMF','lbMSELR','ubMSELR', ...
            'lbMSEMF','ubMSEMF');
    end

end

% matlabpool close % parallel computing

end