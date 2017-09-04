function [] = getCoKrigVar(pah)
% GETCOKRIGVAR This function will calculate and display covariance parameters for PM,
% PAH, and PM/PAH
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

% default variables
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

% reformatting data of PM and PAH
[ZhPM,cMSPM,tMEPM,nanratio]=valstv2stg(chPM,zPM);
[ZhPAH,cMSPAH,tMEPAH,nanratio]=valstv2stg(chjustPAH,zjustPAH);

% in order to come up with rLag, I need to know what the spatial
% distances are and to equal spacing in a log space
cMS = cMSPAH;
tME = tMEPAH;
DMS = sqrt(bsxfun(@plus,dot(cMS,cMS,2),dot(cMS,cMS,2)')-2*(cMS*cMS'));

% spatial lags and tolerances
rLag = prctile(unique(DMS(:)),[0:5:50]);
rTol = [0 (rLag(2:end)-rLag(1:end-1))/2];

% temporal lags and tolerances
ttemp = unique(tME(:))-min(tME(:));
tLag = prctile(ttemp,0:5:50);
tTol = [0 (tLag(2:end)-tLag(1:end-1))/2];

% covariance for PM
[CrPM nprPM]=stcov(ZhPM,cMSPM,tMEPM,ZhPM,cMSPM,tMEPM,rLag,rTol,0,0);
[CtPM nptPM]=stcov(ZhPM,cMSPM,tMEPM,ZhPM,cMSPM,tMEPM,0,0,tLag,tTol);

% covariance for PAH
[CrPAH nprPAH]=stcov(ZhPAH,cMSPAH,tMEPAH,ZhPAH,cMSPAH,tMEPAH,rLag,rTol,0,0);
[CtPAH nptPAH]=stcov(ZhPAH,cMSPAH,tMEPAH,ZhPAH,cMSPAH,tMEPAH,0,0,tLag,tTol);

% cross covariance between PM and PAH
[CrX nprX]=stcov(ZhPM,cMSPM,tMEPM,ZhPAH,cMSPAH,tMEPAH,rLag,rTol,0,0);
[CtX nptX]=stcov(ZhPM,cMSPM,tMEPM,ZhPAH,cMSPAH,tMEPAH,0,0,tLag,tTol);

% spatial: find covariance parameter
f = @(ar1) fitfun(CrPM,CrPAH,CrX,ar1,CrPM(1),CrPAH(1),CrX(1),rLag');  
[arf,fval,exitflag,output] = fminsearch(f,750000);

% spatial: find covariance parameter
f = @(at1) fitfun(CtPM,CtPAH,CtX,at1,CtPM(1),CtPAH(1),CtX(1),tLag);  
[atf,fval,exitflag,output] = fminsearch(f,300);

% covariance model
Cr = [ CrPAH(1) CrX(1) ; CrX(1) CrPM(1) ];
covmodel={'exponentialC/exponentialC'};
covparam={ { [Cr] [arf atf] } };

% saving data
save(sprintf('matfiles/covmodel_cokrig%s.mat',valname{pah+4}),'Cr','covmodel','covparam');

% plotting the spatial range found
figure; hold on;
r=0:1000:300000;
subplot(2,2,1);
plot(rLag,CrPM,'bo');
title('Covariance PM');
hold on;
plot(r,CrPM(1).*exp(-3.*r./arf),'r-');
ylabel('C(r,t=0 days)');
xlabel(sprintf('Spatial lag r (meters), ar = %0.2f',arf));

subplot(2,2,2); 
plot(rLag,CrX,'bo');
title(sprintf('Cross-Covariance PM-PAH: %s',valdispname{pah}));
hold on;
plot(r,CrX(1).*exp(-3.*r./arf),'r-');
ylabel('C(r,t=0 days)');
xlabel(sprintf('Spatial lag r (meters), ar = %0.2f',arf));

subplot(2,2,4); hold on;
plot(rLag,CrPAH,'bo');
title(sprintf('Covariance PAH: %s',valdispname{pah}));
plot(r,CrPAH(1).*exp(-3.*r./arf),'r-');
ylabel('C(r,t=0 days)');
xlabel(sprintf('Spatial lag r (meters), ar = %0.2f',arf));

% save figure 
set(gcf,'Position',[0 0 800 500]); 
print(gcf,'-painters','-dpng','-r600',sprintf('figures/covmodel_cokrig_space_%s.png',valname{pah+4}));

% plotting the temporal range found
figure; hold on;
t=0:1:250;
subplot(2,2,1);
plot(tLag,CtPM,'bo');
title('Covariance PM');
hold on;
plot(t,CtPM(1).*exp(-3.*t./atf),'r-');
ylabel('C(r=0 meters,t)');
xlabel(sprintf('Temporal lag t (days), at = %0.2f',atf));

subplot(2,2,2); 
plot(tLag,CtX,'bo');
title(sprintf('Cross-Covariance PM-PAH: %s',valdispname{pah}));
hold on;
plot(t,CtX(1).*exp(-3.*t./atf),'r-');
ylabel('C(r=0 meters,t)');
xlabel(sprintf('Temporal lag t (days), at = %0.2f',atf));

subplot(2,2,4); hold on;
plot(tLag,CtPAH,'bo');
title(sprintf('Covariance PAH: %s',valdispname{pah}));
plot(t,CtPAH(1).*exp(-3.*t./atf),'r-');
ylabel('C(r=0 meters,t)');
xlabel(sprintf('Temporal lag t (days), at = %0.2f',atf));

% save figure 
set(gcf,'Position',[0 0 800 500]); 
print(gcf,'-painters','-dpng','-r600',sprintf('figures/covmodel_cokrig_time_%s.png',valname{pah+4}));

end