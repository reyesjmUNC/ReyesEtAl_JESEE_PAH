function [] = autoCovModel_PM()
% AUTOCOVMODEL_PM automating the covariance model for PM2.5

% load data
load('matfiles/pah_data.mat');
zPM = val(:,4);
cPAH = [ProjectX ProjectY Time];

% reformatting data
[Zh,cMS,tME,nanratio]=valstv2stg(cPAH,zPM);

% in order to come up with rLag, I need to know what the spatial
% distances are and to equal spacing in a log space
DMS = sqrt(bsxfun(@plus,dot(cMS,cMS,2),dot(cMS,cMS,2)')-2*(cMS*cMS'));

% spatial lags and tolerances
rLag = [0 prctile(unique(DMS(:)),[2:10 12.5 15:5:50])];
rTol = [0 (rLag(2:end)-rLag(1:end-1))/2];

% temporal lags and tolerances
ttemp = unique(tME(:))-min(tME(:));
tLag = [0 prctile(unique(ttemp(:)),[0.25 0.5 0.75 1 1.5 2:10 12.5 15])];
tTol = [0 (tLag(2:end)-tLag(1:end-1))/2];

% find experimental covariance
[Cr npr]=stcov(Zh,cMS,tME,Zh,cMS,tME,rLag,rTol,0,0);
[Ct npt]=stcov(Zh,cMS,tME,Zh,cMS,tME,0,0,tLag,tTol);
idxr = ~isnan(Cr) & Cr>0;
Crtest = Cr(idxr); rLagtest = rLag(idxr);
idxt = ~isnan(Ct) & Ct>0;
Cttest = Ct(idxt); tLagtest = tLag(idxt);
Cr1 = Cr(1);

% finding joint space/time covariance model
% joint exponential exponential 
s = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0,0,0,0],...
    'Upper',[1,Inf,Inf,Inf,Inf],'Startpoint',[0.5,200000,2000000,500,300]);
g = fittype( 'jointexpexp(alp,ar1,ar2,at1,at2,Cr1,x,y)','options',s,...
    'problem',{'Cr1'},'independent',{'x','y'},'dependent',{'z'} );

x = [rLag' ; zeros(length(tLag),1)];
y = [zeros(length(rLag),1) ; tLag'];
z = [Cr ; Ct'];
[f gof output] = fit([x,y],z,g,'problem',Cr1);
disp(f);

covmodel = {'exponentialC/exponentialC','exponentialC/exponentialC'};
covparam = {[f.alp*Cr(1) f.ar1 f.at1],[(1-f.alp)*Cr(1) f.ar2 f.at2]};
save('matfiles/covmodel_PM2p5.mat','covmodel','covparam', ...
    'f','gof','rLag','output','s','g','tLag','Ct','npt','Cr','npr');

% plotting the spatial covariance
r=0:1000:rLag(end);
subplot(2,1,1);
plot(rLag,Cr,'bo');
hold on;
plot(r,Cr(1).*(f.alp.*exp(-3.*r./f.ar1) + (1-f.alp).*exp(-3*r./f.ar2)),'r-');
ylabel('C(r,t=0 days)');
xlabel(sprintf('spatial lag r (meters), ar1 = %f, ar2 = %f',f.ar1,f.ar2));
title('Covariance model for PM2.5');

% plotting the temporal covariance
t=0:1:250;
subplot(2,1,2);
plot(tLag,Ct,'bo');
hold on;
plot(t,Ct(1).*(f.alp.*exp(-3.*t./f.at1) + (1-f.alp).*exp(-3*t./f.at2)),'r-');
ylabel('C(r=0 meters,t)');
xlabel(sprintf('temporal lag t (days), at1 = %f, at2 = %f',f.at1,f.at2));

% save figure 
set(gcf,'Position',[0 0 800 500]); 
print(gcf,'-painters','-dpng','-r600','figures/covmodel_PM2p5.png');

end