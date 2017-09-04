function [Z,cMS,tME,nanratio]=valstv2stg_new(p,z)
% valstv2stg                - Converts values from s/t vector coord to s/t grid coord 
%
% Converts the values of a space/time variable from a s/t vector 
% format (i.e. the variable z is listed as a vector of n values)
% to a grid format (i.e. the variable Z is given as a nMS by nME matrix 
% corresponding to nMS Monitoring Sites and nME Measuring Events).
% Use help stgridsyntax for information on the s/t grid format.
%
% SYNTAX :
%
% [Z,cMS,tME,nanratio]=valstv2stg(ch,z,cMS,tME);
%
% INPUT :
%
% ch         n by d+1     matrix of space/time coordinates for spatial domain of dimension d
% z          n by 1       vector of field value at coordinate ch
%
% OUTPUT :
%
% Z          nMS by nME   matrix of values for the variable Z corresponding to 
%                         nMS Monitoring Sites and nME Measuring Event
% cMS        nMS by d     matrix of spatial coordinates for the nMS Measuring Sites
% tME        1 by nME     vector of times of the tME Measuring Events
% nanratio   scalar       ratio of the NaNs in Z (0<=nanratio<=1) 
%
% 
% See also stgridsyntax for help on s/t grid format
 
[n,m]=size(p);
if m<2, error('p must have at least two columns'); end
d=m-1;
if size(z,1)~=n, error('p and z must have the same number of rows'); end
if size(z,2)~=1, error('z must have only one column'); end

if isdupli(p)
disp('warning in valstv2stg: duplicated data were averaged');
[p,z] = avedupli(p,z);
end

cMS=unique(p(:,1:d),'rows');
nMS=size(cMS,1);
tME=unique(p(:,end))';
nME=length(tME);  

Z(1:nMS,1:nME)=NaN;
s = [nMS,nME];
[i j] =ind2sub(s,1:nMS*nME);
Zv = [cMS(i,:) tME(j)'];
[C ia ib] = intersect(Zv,p,'rows');
Z(ia) = z(ib);

nanratio=sum(isnan(Z(:)))/(nMS*nME);

end
