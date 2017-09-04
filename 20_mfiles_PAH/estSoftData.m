function [] = estSoftData(pah,LR)
%
% ESTSOFTDATA Determines Mean and Variance of Soft PAH data
% This function determines the mean and variances for PAH at specified
% nPAHs and dmax values determined by the PAH name in the input.
%
% inputs:   pah      - This is the name of the PAH being analyzed. This file
%                      should be in single quotes and have the same name as
%                      listed in the heading of PAHData.txt, ie 'chrysene'.
%                      Default: 'benz_a_anthracene'
%                      Type: string
%           LR       - This variable determines which method is used to
%                      find the soft data. If LR = 1, a linear regression 
%                      method will be used to determine the soft data. If 
%                      LR = 0, a mass fraction apprach will be used to 
%                      determine the soft data.
%                      Default: 1
%                      Size: 1x1
%
% outputs:  
% There are no outputs.
%

% default variables
if nargin < 1, pah = 1; end
if nargin < 2, LR = 0; end

% defining strings
if LR == 1, LRStr = 'LR'; else LRStr = 'MF'; end

% load data
load('matfiles/pah_data.mat');
idxPAH = ~isnan(val(:,pah+4));
load(sprintf('matfiles/bestMinNeigh_%s_ellip_%s.mat',valname{pah+4},LRStr));

% all the raw data we need
cPM = [ProjectX ProjectY Time];
cPAH = [ProjectX(idxPAH) ProjectY(idxPAH) Time(idxPAH)];
zPM = log(val(:,4));
zPMsub = log(val(idxPAH,4));
zPAH = log(val(idxPAH,pah+4));
zPAHall = log(val(:,pah+4)); zPAHall(~idxPAH) = NaN;
zMF = log(val(idxPAH,pah+4)./val(idxPAH,4));

% calculating mean and variance
mPAHs = NaN*ones(length(zPAHall),1);
vPAHs = NaN*ones(length(zPAHall),1);
csub = cell(length(zPAHall),1);
Zsub = cell(length(zPAHall),1);
dsub = cell(length(zPAHall),1);
nsub = cell(length(zPAHall),1);
index = cell(length(zPAHall),1);
PDFpnts = cell(length(zPAHall),1);
            
for i = 1:length(idxPAH)

    if isnan(zPAHall(i))
        if LR == 1 % linear regression        
            
            % for ellipsoidal neighborhood        
            nPAHs = minR;
            dmax = [ 1000000 1000 minC ];

            [csub{i},Zsub{i},dsub{i},nsub{i},index{i}]=neighbours(cPM(i,:),cPAH,zPAH,nPAHs,dmax);
            if length(index{i}) > 1
                LRs = regstats(Zsub{i},zPMsub(index{i}),'linear');
                mPAHs(i) = LRs.beta(1)+LRs.beta(2)*zPM(i);
                % calculate variance (prediction interval)
                [p,S,mu] = polyfit(zPMsub(index{i}),Zsub{i},1);
                [Yhat,DELTA] = polyconf(p,zPM(i),S);
                if DELTA == Inf, DELTA = 1000; end % pick a big number here
                vPAHs(i) = (DELTA./2).^2; 
                PDFpnts{i} = LRs.r+LRs.beta(1)+LRs.beta(2)*zPM(i);
            elseif length(index{i}) == 1
                mPAHs(i) = mean(Zsub{i})+zPM(i);
                vPAHs(i) = var(zPAH);
                PDFpnts{i} = Zsub{i}+zPM(i);
                % if only one data point in the cylinder, assign that
                % PAH value and pick an arbitrary variance
            elseif length(index{i}) == 0
                [csub{i},Zsub{i},dsub{i},nsub{i},index{i}]=neighbours(cPM(i,:),cPAH,zPAH,1,[1000000 1000 1]);
                mPAHs(i) = mean(Zsub{i}) + zPM(i);
                vPAHs(i) = var(zPAH);
                PDFpnts{i} = Zsub{i}+zPM(i);
                % if no data is in the cylinder, pick the closest point
                % and pick an arbitrary variance
            end            
            
        elseif LR == 0 % mass fraction

            % for ellipsoidal neighborhood        
            nPAHs = minR;
            dmax = [ 1000000 1000 minC ];

            [csub{i},Zsub{i},dsub{i},nsub{i},index{i}]=neighbours(cPM(i,:),cPAH,zMF,nPAHs,dmax);  
            if length(index{i}) > 0
                mPAHs(i) = mean(Zsub{i})+zPM(i);
                vPAHs(i) = var(Zsub{i}); 
                PDFpnts{i} = Zsub{i}+zPM(i);
            elseif length(index{i}) == 0
                [csub{i},Zsub{i},dsub{i},nsub{i},index{i}]=neighbours(cPM(i,:),cPAH,zMF,1,[1000000 1000 1]);
                mPAHs(i) = mean(Zsub{i})+zPM(i);
                vPAHs(i) = var(zPAH); 
                PDFpnts{i} = Zsub{i}+zPM(i);
                % if no data is in the cylinder, pick the closest point
                % and pick an arbitrary variance
            end
            
        end            
    else
        
        mPAHs(i) = zPAHall(i);
        vPAHs(i) = 0;
        
    end
    
end

% saving soft data
save(sprintf('matfiles/soft_%s_ellip_%s.mat',valname{pah+4},LRStr), ...
    'cPM','mPAHs','vPAHs','csub','Zsub','dsub','nsub','index','PDFpnts', ...
    'cPAH','zPAH','zPM','zPMsub','zPAHall','zMF'); 

end