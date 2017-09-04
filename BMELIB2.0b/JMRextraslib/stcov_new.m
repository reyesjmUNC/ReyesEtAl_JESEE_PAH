function [C,np]=stcov_new(Zi,cMSi,tMEi,Zj,cMSj,tMEj,rLag,rLagTol,tLag,tLagTol)

% stcov                     - Estimate space/time cross covariance values from data
% 
% Estimates the space/time cross covariance between variable Zi and
% variable Zj with measurements at fixed monitoring stations. The 
% monitoring stations for Zi are located at coordinates cMSi, and 
% measuring events are at times tMEi, and similarly cMSj and tMEj are
% the coordinates of monitoring sites and times of measuring events
% for Zj.
%
% SYNTAX :
%
% [C np]=stcov(Zi,cMSi,tMEi,Zj,cMSj,tMEj,rLag,rLagTol,tLag,tLagTol)
%
% INPUT :
%
% Zi      nMSi by nMEi matrix with measurements of Zi at the nMSi monitoring
%                      sites and nMEi measuring events. Zi may have NaN values.
% cMSi    nMSi by 2    matrix of spatial x-y coordinates for the nMSi monitoring
%                      sites
% tMEi    1 by nMEi    vector with the time of the measuring events
% Zj      nMSj by nMEj matrix with measurements of Zj at the nMSj monitoring
%                      sites and nMEj measuring events. Zj may have NaN values.
% cMSj    nMSj by 2    matrix of spatial x-y coordinates for the nMSj monitoring
%                      sites
% tMEj    1 by nMEj    vector with the time of the measuring events
% rLag    nr by 1      vector with the r lags
% rLagTol nr by 1      vector with the tolerance for the r lags
% tLag    1 by nt      vector with the t lags
% tLagTol nt by 1      vector with the tolerance for the t lags
%
% OUTPUT :
%
% C      nr by nt      matrix with the covariance values at lags rLag
%                      and tlags
% np     nr by nt      matrix with the number of pairs
%
% NOTE : 
% 
% Use help stgridsyntax for help on s/t grid format

%
% Verify input
%
if size(cMSi,1)~=size(Zi,1) || size(cMSi,2)~=2
    error('cMSi must be a nMSi by 2 matrix'); 
end;
if size(tMEi,1)~=1 || size(tMEi,2)~=size(Zi,2)
  error('tMEi must be a 1 by nMEi vector'); 
end;
if size(cMSj,1)~=size(Zj,1) || size(cMSj,2)~=2
  error('cMSj must be a nMSj by 2 matrix'); 
end;
if size(tMEj,1)~=1 || size(tMEj,2)~=size(Zj,2)
  error('tMEj must be a 1 by nMEj vector'); 
end;
nr=length(rLag);
if length(rLagTol)~=nr
  error('rLag and rLagTol must have same length');
end;
nt=length(tLag);
if length(tLagTol)~=nt
  error('tLag and tLagTol must have same length');
end;

% all the distances in space 
DMS = sqrt(bsxfun(@plus,dot(cMSi,cMSi,2),dot(cMSj,cMSj,2)')-2*(cMSi*cMSj'));
rowr = cell(nr,1);
colr = cell(nr,1);
for ir = 1:nr
    [rowr{ir} colr{ir}] = find( DMS>=rLag(ir)-rLagTol(ir) & DMS<=rLag(ir)+rLagTol(ir) );
end

% all the distances in time
DME = abs(bsxfun(@minus,tMEi,tMEj'))';
rowt = cell(nt,1);
colt = cell(nt,1);
for it = 1:nt
    [rowt{it} colt{it}] = find( DME>=tLag(it)-tLagTol(it) & DME<=tLag(it)+tLagTol(it) ); 
end
               
% covariance calculation
irs = repmat([1:nr]',nt,1);
its = repmat([1:nt],nr,1); its = its(:);
irt = [irs its];

C = NaN*ones(nr*nt,1);
np = zeros(nr*nt,1);
nrnt = nr*nt;

for i = 1:nrnt
    
    Zr = Zi(rowr{irt(i,1)},rowt{irt(i,2)}); 
    Zr = Zr(:);
    Zc = Zj(colr{irt(i,1)},colt{irt(i,2)}); 
    Zc = Zc(:);
    idx = ~isnan(Zr) & ~isnan(Zc);
    sumidx = sum(idx);

    if sumidx>0, 
        Zridx = Zr(idx);
        Zcidx = Zc(idx);
        C(i) = mean2(Zridx.*Zcidx) - mean2(Zridx)*mean2(Zcidx);
        np(i) = sumidx; 
    end 
    
end

C = reshape(C,nr,nt);
np = reshape(np,nr,nt);

end
