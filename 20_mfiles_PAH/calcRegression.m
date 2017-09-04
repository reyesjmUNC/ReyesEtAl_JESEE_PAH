function [MSELR,AICLR,MSEMF,AICMF,lbMSELR,ubMSELR,lbMSEMF,ubMSEMF] = ...
    calcRegression(ch,PM,PAH,spatrad,temprad,nmax,stmet)

if nargin < 1, ch = rand(10,3); end
if nargin < 2, PM = rand(10,1); end
if nargin < 3, PAH = rand(10,1); end
if nargin < 4, spatrad = 100; end
if nargin < 5, temprad = 100; end
if nargin < 6, nmax = 90; end
if nargin < 7, stmet = 0.5; end

idxPAH = ~isnan(PAH);
ch = ch(idxPAH,:); PM = log(PM(idxPAH)); PAH = log(PAH(idxPAH));
predLRs = NaN*ones(length(ch),1);
AICLRs = NaN*ones(length(ch),1);
predMFs = NaN*ones(length(ch),1);
AICMFs = NaN*ones(length(ch),1);

for i = 1:length(ch)
    
    % calculating all space/time distances from the point
    ds = sqrt( (ch(:,1)-ch(i,1)).^2 + (ch(:,2)-ch(i,2)).^2 );
    dt = abs( ch(:,3)-ch(i,3) );   
    idx = ds <= spatrad & dt <= temprad;
    idx(i) = 0; % takes out the point being evaluated
    num = sum(idx); % number of points in the space/time cylinder
    
    % adding the nmax condition
    if num > nmax
        dmax = [ spatrad temprad stmet ];
        [csub,Zsub,dsub,nsub,idx] = neighbours2(ch(i,:),ch(idx,:),PAH(idx),nmax,dmax); 
        temp = zeros(length(ch),1); temp(idx) = 1; idx = logical(temp);
        num = nmax;
    end
      
    % getting the MSE and AIC value for the LR model
    if sum(idx) > 1 % if there are points
        stats = regstats(PAH(idx),PM(idx));
        predLRs(i) = PM(i)*stats.beta(2) + stats.beta(1);
        AICLRs(i) = 2 + num*log( (predLRs(i) - PAH(i)).^2./num );

        % getting the MSE and AIC value for the MF model
        meanMF = mean(PAH(idx)-PM(idx));
        predMFs(i) = PM(i) + meanMF;
        AICMFs(i) = 4 + num*log( (predMFs(i) - PAH(i)).^2./num );
    end
    
end

n = length((predLRs-PAH).^2);
MSELR = mean((predLRs-PAH).^2);
lbMSELR = mean((predLRs-PAH).^2) - 2*sqrt( var((predLRs-PAH).^2)/n );
ubMSELR = mean((predLRs-PAH).^2) + 2*sqrt( var((predLRs-PAH).^2)/n );
AICLR = mean(AICLRs);
MSEMF = mean((predMFs-PAH).^2);
lbMSEMF = mean((predMFs-PAH).^2) - 2*sqrt( var((predMFs-PAH).^2)/n );
ubMSEMF = mean((predMFs-PAH).^2) + 2*sqrt( var((predMFs-PAH).^2)/n );
AICMF = mean(AICMFs);

end