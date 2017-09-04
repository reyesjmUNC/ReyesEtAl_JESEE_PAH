function [] = autoCovmodel(pah,mf)
% AUTOCOVMODEL automating the covariance model
%
% INPUT:
%
% pahname   k by 1      a string of the name of the pah being evaluated
% mf        1 by 1      a binary variable that states if the mass fraction is being
%                       calculated If yes, mf = 1.
%               
% OUTPUT: There are no outputs.
%

if nargin < 1, pah = 1; end
if nargin < 2, mf = 0; end
if mf == 1 , 
    mfStr = '_mf'; 
    mfStrdisp = 'mass fraction';
else
    mfStr = '';
    mfStrdisp = '';
end

% load data
load('matfiles/pah_data.mat');
idx = ~isnan(val(:,pah+4));
zPAH = val(idx,pah+4);
zPM = val(idx,4);
cPAH = [ProjectX(idx) ProjectY(idx) Time(idx)];

% determining mass fraction
if mf == 1
    zPAH = log(zPAH./zPM);
else
    zPAH = log(zPAH);
end
    
% reformatting data
[Zh,cMS,tME,nanratio]=valstv2stg(cPAH,zPAH);

% in order to come up with rLag, I need to know what the spatial
% distances are and to equal spacing in a log space
DMS = sqrt(bsxfun(@plus,dot(cMS,cMS,2),dot(cMS,cMS,2)')-2*(cMS*cMS'));

% spatial lags and tolerances
rLag = prctile(unique(DMS(:)),[0:5:50]);
rTol = [0 (rLag(2:end)-rLag(1:end-1))/2];

% temporal lags and tolerances
ttemp = unique(tME(:))-min(tME(:));
tLag = prctile(ttemp,0:5:50);
tTol = [0 (tLag(2:end)-tLag(1:end-1))/2];

% find experimental covariance
[Cr npr]=stcov(Zh,cMS,tME,Zh,cMS,tME,rLag,rTol,0,0);
[Ct npt]=stcov(Zh,cMS,tME,Zh,cMS,tME,0,0,tLag,tTol);

%%% spatial: find covariance model/parameters
s = fitoptions('Method','NonlinearLeastSquares','Lower',[0], ...
    'Upper',[Inf],'Startpoint',[750000]);
g = fittype( 'expmodel(ar1,Cr1,x)','options',s,...
    'problem',{'Cr1'},'independent',{'x'},'dependent',{'z'} );

idxnan = ~isnan(Cr) & Cr>0;
x = rLag(idxnan);
z = Cr(idxnan);
Cr1 = z(1);
[f gof output] = fit(x',z,g,'problem',Cr1); 
arf = f.ar1;
disp(f);

%%% temporal: find covariance model/parameters
s = fitoptions('Method','NonlinearLeastSquares','Lower',[0], ...
    'Upper',[Inf],'Startpoint',[60]);
g = fittype( 'expmodel(ar1,Cr1,x)','options',s,...
    'problem',{'Cr1'},'independent',{'x'},'dependent',{'z'} );

idxnan = ~isnan(Ct) & Ct>0 ;
x = tLag(idxnan);
z = Ct(idxnan);
Ct1 = z(1);
[f gof output] = fit(x',z',g,'problem',Ct1); 
atf = f.ar1;
disp(f);

% covariance model
covmodel={'exponentialC/exponentialC'};
covparam={[Cr(1) arf atf]};

% saving data
save(sprintf('matfiles/covmodel_%s%s.mat',valname{pah+4},mfStr),'covmodel','covparam','tLag', ...
    'Ct','npt','rLag','Cr','npr','arf','atf');

figure; hold on;

% plotting the spatial covariance
r=0:1000:300000;
subplot(2,1,1);
plot(rLag,Cr,'bo');
hold on;
plot(r,Cr(1).*exp(-3.*r./arf),'r-');
ylabel('C(r,t=0 days)');
xlabel(sprintf('Spatial lag r (meters), ar = %f',arf));
title(sprintf('Covariance model for %s %s',valdispname{pah},mfStrdisp));

% plotting the temporal covariance
t=0:1:250;
subplot(2,1,2);
plot(tLag,Ct,'bo');
hold on;
plot(t,Ct(1).*exp(-3.*t./atf),'r-');
ylabel('C(r=0 meters,t)');
xlabel(sprintf('temporal lag t (days), at = %f',atf));

% save figure 
set(gcf,'Position',[0 0 800 500]); 
print(gcf,'-painters','-dpng','-r600',sprintf('figures/covmodel_%s%s.png',valname{pah+4},mfStr));

end
