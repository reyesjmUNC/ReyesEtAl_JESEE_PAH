function [zk,vk,temp1,temp2,temp1a,temp1b,allhard,alllambda1,alllambda2,allch,allcs]= ...
    krigingME2_correct_parallel(ck,ch,cs,zh,zs,vs,covmodel,covparam,nhmax,nsmax,dmax,order,options);

% krigingME                 - prediction using kriging with measurement errors (Jan 1,2001)
%
% Account for uncertainties by using a variation of kriging.
% The uncertainties at the data locations are modeled as random
% noises, and the resulting data can be considered as probabilitic
% soft data. These probabilistic soft data are assumed to be
% completely characterized by their mean and their variance. As
% these are only second-order moments related information, a linear
% kriging algorithm can still process adequately any mixture of these
% soft data with hard data. 
% 
% SYNTAX :
%
% [zk,vk]=krigingME(ck,ch,cs,zh,zs,vs,covmodel,covparam,nhmax,nsmax,dmax,order,options);
%
% INPUT :
%
% ck         nk by d   matrix of coordinates for the estimation locations.
%                      A line corresponds to the vector of coordinates at
%                      an estimation location, so the number of columns
%                      corresponds to the dimension of the space. There is
%                      no restriction on the dimension of the space.
% ch         nh by d   matrix of coordinates for the hard data locations,
%                      with the same convention as for ck.
% cs         ns by d   matrix of coordinates for the soft data locations,
%                      with the same convention as for ck.
% zh         nh by 1   vector of values for the hard data at the coordinates
%                      specified in ch.
% zs         ns by 1   vector of values for the mean of the soft data at the
%                      coordinates specified in cs.
% vs         ns by 1   vector of values for the variance of the soft data at
%                      the coordinates specified in cs.
% covmodel   string    string that contains the name of the covariance model
%                      that is used for the estimation (see the MODELS directory).
%                      Variogram models can not be used with this function. 
% covparam   1 by k    vector of values for the parameters of covmodel, according
%                      to the convention for the corresponding covariance model. 
% nhmax      scalar    maximum number of hard data values that are considered
%                      for the estimation at the locations specified in ck.
% nsmax      scalar    maximum number of soft data values that are considered
%                      for the estimation at the locations specified in ck.
% dmax       scalar    maximum distance between an estimation location and
%                      existing hard/soft data locations. All hard/soft data
%                      locations separated by a distance smaller than dmax from
%                      an estimation location will be included in the estimation
%                      process for that location, whereas other hard/soft data
%                      locations are neglected.
% order      scalar    order of the polynomial mean along the spatial axes at the
%                      estimation locations. For the zero-mean case, NaN (Not-a-
%                      Number) is used.
% options    scalar    optional parameter that can be used if the default value
%                      is not satisfactory (otherwise it can simply be omitted
%                      from the input list of variables). options(1) is taking
%                      the value 1 or 0 if the user wants or does not want to
%                      display the order number of the location which is
%                      currently processed, respectively.
%
% OUTPUT :
%
% zk         nk by 1   vector of estimated values at the estimation locations. A
%                      value coded as NaN means that no estimation has been performed
%                      at that location due to the lack of available data. 
% vk         nk by 1   vector of estimation (kriging) variances at the estimation
%                      locations. As for zk, a value coded as NaN means that no
%                      estimation has been performed at the corresponding location.
%
% NOTE :
%
% 1- Note that in the case there are no available hard data at all,
% ch and zh can be entered as the empty [ ] matrices.
%
% 2- All the specific conventions for specifying nested models,
% multivariate or space-time cases are the same as for kriging.m.

%%%%%% Initialize the parameters

if nargin<13,
  options(1)=0;
end;

noindex=~iscell(ck);       % test if there is an index for the variables
if noindex==1,
  nk=size(ck,1);           % nk is the number of estimation points 
  nh=size(ch,1);           % nh is the number of hard data
  ns=size(cs,1);           % ns is the number of soft data
else
  nk=size(ck{1},1);
  nindexk=length(ck{2});
  if nindexk==1,
    ck{2}=ck{2}*ones(nk,1);
  end;
  nh=size(ch{1},1);
  ns=size(cs{1},1);
end;

if options(1)==1,
  num2strnk=num2str(nk);
end;

zk=zeros(nk,1)*NaN;
vk=zeros(nk,1)*NaN;
temp1=zeros(nk,1)*NaN;
temp2=zeros(nk,1)*NaN;
temp1a=zeros(nk,1)*NaN;
temp1b=zeros(nk,1)*NaN;
allhard = cell(nk,1);
alllambda1 = cell(nk,1);
alllambda2 = cell(nk,1);
allch = cell(nk,1);
allcs = cell(nk,1);

%%%%%% Main loop starts here
%matlabpool open 12
parfor i=1:nk
%for i=1:nk   
  if mod(i,1000) == 0, disp(i); end
  %disp(i);
  if noindex==1,
    ck0=ck(i,:);
  else
    ck0={ck{1}(i,:),ck{2}(i)};
  end;

  [chlocal,zhlocal,dh,sumnhlocal,index]=neighbours2a(ck0,ch,zh,nhmax,dmax);
  [cslocal,zslocal,dh,sumnslocal,index]=neighbours2a(ck0,cs,zs,nsmax,dmax);
  vslocal=vs(index);
  
  allhard{i} = zhlocal;
  alllambda1{i} = zslocal;
  alllambda2{i} = vslocal;
  allch{i} = chlocal;
  allcs{i} = cslocal;

  Khh=coord2K(chlocal,chlocal,covmodel,covparam);      % built the left-hand side matrix for hard data
  Kss=coord2K(cslocal,cslocal,covmodel,covparam);      % built the left-hand side matrix for soft data
  Khs=coord2K(chlocal,cslocal,covmodel,covparam);      % built the left-hand side matrix for hard-soft data
  Kss=Kss+diag(vslocal);                               % add the error variances on the diagonal
  kh=coord2K(chlocal,ck0,covmodel,covparam);           % built the right-hand side vector for hard data
  ks=coord2K(cslocal,ck0,covmodel,covparam);           % built the right-hand side vector for soft data

  K=[[Khh,Khs];[Khs',Kss]];                            % built the composite left-hand side matrix
  k=[kh;ks];                                           % built the composite right-hand side matrix

  if noindex==1,
    chslocal=[chlocal;cslocal];
  else
    chslocal{1}=[chlocal{1};cslocal{1}];
    chslocal{2}=[chlocal{2};cslocal{2}];
  end;
  [X,x]=krigconstr(chslocal,ck0,order);                % build the constraint matrices

  index=findpairs(ck0,cslocal);                        % test if there is a soft data at estimation point
  %if isempty(index),
    k0=coord2K(ck0,ck0,covmodel,covparam);             % compute the variance at ck0
  %else
  %  k0=coord2K(ck0,ck0,covmodel,covparam)+vslocal(index(2));
  %end;
  if (sumnhlocal+sumnslocal)>0,                        % if there is at least one hard or soft data
    nx=size(X,2);
    Kadd=[[K,X];[X',zeros(nx)]];
    kadd=[k;x];
    lam=Kadd\kadd;                                     % compute the kriging weights lam
    lam=lam(1:(sumnhlocal+sumnslocal));                % remove the Lagrangians from the solution
    lamt=lam';
    
    avgnegativeweight = mean(abs(lamt(lamt<0)));
    avgcov = mean(kadd(lam<0));
    weightformean = 1-sum(lamt);
    
    newlamt = lamt;
    newlamt(lamt < 0) = 0;
    newlamt((lamt > 0) & (kadd' < avgcov) & (lamt < avgnegativeweight)) = 0;
    newweight = [newlamt,weightformean]./sum([newlamt,weightformean]);
    newweight = newweight(1:end-1);
    zk(i)=newweight*[zhlocal;zslocal];                      % compute the kriging estimates zk(i)
    vk(i)=k0-2*newweight*k+newweight*K*newweight';          % compute the kriging variance vk(i)
    newweighth = newweight(1:sumnhlocal); newweights = newweight(sumnhlocal+1:end);
    
    if sum(newweighth==0)==length(newweighth)  
        temp1(i) = NaN;
        temp1a(i) = NaN;
        temp1b(i) = NaN;
    else
        temp1(i) = nanmean(zhlocal(newweighth~=0));
        temp1a(i) = nanmax(zhlocal(newweighth~=0));
        temp1b(i) = nanmin(zhlocal(newweighth~=0));
    end
    if sum(newweights==0)==length(newweights) 
        temp2(i) = NaN;
    else
        temp2(i) = nanmean(zslocal(newweights~=0));
    end
    
  end;

%   if options(1)==1,
%     disp([num2str(i),'/',num2strnk]);
%   end;
end;
%matlabpool close

