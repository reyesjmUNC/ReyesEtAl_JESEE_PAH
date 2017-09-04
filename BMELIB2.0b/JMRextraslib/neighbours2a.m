function [csub,Zsub,dsub,nsub,index]=neighbours2a(c0,c,Z,nmax,dmax);% neighbours                - radial neighbourhood selection (Jan 1,2001)%% Select a subset of coordinates and variables based on% their distances from the coordinate c0.%% SYNTAX :%% [csub,Zsub,dsub,nsub,index]=neighbours(c0,c,Z,nmax,dmax);%% INPUT :%% c0      1 by d   vector of coordinates, where d is the dimension%                  of the space% c       n by d   matrix of coordinates% Z       n by k   matrix of values, where each column corresponds to%                  the values of a same variable and each line corresponds%                  to the values of the diferent variables at the corresponding%                  c coordinates% nmax    scalar     maximum number of lines of Z that must be kept.% dmax    scalar     maximum Euclidian distance between c0 and c.%% OUTPUT :%% csub    m by d   matrix which is a subset of lines of the c matrix (m<=n)% Zsub    m by k   matrix which is a subset of lines of the Z matrix.% dsub    m by 1   vector of distances between csub and c0.% nsub    scalar   length of the dsub vector.% index   m by 1   vector giving the ordering of the lines in csub with%                  respect to the initial matrix c.%% NOTE :%% 1-In the case of space/time coordinates, dmax is a vector of length 3,% and the last column of c0 and c is the temporal coordinate. In this case,% dmax(1) is the maximum spatial distance between the coordinate in c and c0,% dmax(2) is the maximum temporal distance between coordinate in c and c0,% and dmax(3) refers to the space/time metric, such that :%% space/time distance=spatial distance+dmax(3)*temporal distance.% % The space/time distance is used to select the nmax closest coordinates c% from the estimation coordinates c0.%% 2- It is possible to specify additional 1 by 1 and n by 1 index vectors,% taking integer values from 1 to nv. The values in the index vectors specify% which of the nv variable is known at each one of the corresponding coordinates.% The c0 and c matrices of coordinates and the index vectors are then grouped% together using the MATLAB cell array notation, so that c0={c0,index0} and% c={c,index}.%%%%%% Test if c is non empty and is cell arrayisST=(length(dmax)==3);if length(dmax)==2,  error('dmax must have 1 element for spatial cases and three elements for space/time cases');end;if iscell(c)==1,  noelements=isempty(c{1});  noindex=0;else   noelements=isempty(c);  noindex=1;end;if noelements==1,  Zsub=[];  dsub=[];  nsub=0;  index=[];  if noindex==1,    csub=[];  else    csub={[],[]};  end;else  %%%%%% two steps data selection inside the neighbourhood   if noindex==1,                   % if there is no index    n=size(c,1);                   % create one with all values=1    c0={c0,1};    c={c,ones(n,1)};  else    n=size(c{1},1);  end;  if ~isST,                          % for the spatial or temporal case,    [d]=coord2dist(c{1},c0{1});      % compute the distances in space or time    index=find(d<=dmax);             % find distances<=dmax  else                                            % for the space-time case    nd=size(c{1},2)-1;                            % compute the dimension of the space    [ds]=sqrt((c{1}(:,1)-c0{1}(:,1)).^2 + (c{1}(:,2)-c0{1}(:,2)).^2); % compute the distances in space    [dt]=abs(c{1}(:,nd+1)-c0{1}(:,nd+1)); % compute the distances in time    index=find((ds<=dmax(1))&(dt<=dmax(2)));      % find distances in space & time<=dmax  end;  %%%%%% checking for the number of data  nv=max(c{2}(index));                % determine maximum value of this index  nsub=0;                             % initialize the value of nsub to 0  indexi=cell(1,nv);  for i=1:nv,                         % loop over the variables    select=find(c{2}(index)==i);      % select coordinates for variable i    indexi{i}=index(select);          % select corresponding values of index    nsub(i)=length(indexi{i});        % compute the number of data for each variable    if nsub(i)>nmax(i),               % when number of data exceeds nmax(i)      if ~isST,        di=d(indexi{i});      else        di=ds(indexi{i})+dmax(3)*dt(indexi{i}); % combine distance in space and time      end;      S=sortrows([indexi{i},di],2);   % sort the distance in increasing order      indexi{i}=S(1:nmax(i),1);       % select the nmax(i) closest data      nsub(i)=nmax(i);    end;  end;  index=[];  for i=1:nv,                         % concatenate the indexi{} as the variable index    index=[index;indexi{i}];          % in increasing order of the index  end;  if ~isST,                           % for the spatial or temporal case,    dsub=d(index);                    % extract the subset of distances in space or time  else                                % for the space-time case,    dsub=[ds(index),dt(index)];       % extract the subset of distances in space and time  end;  nsub=sum(nsub);                     % compute the length of the dsub vector   Zsub=Z(index,:);                    % extract the subset of lines in V   if noindex==1,                      % if there was no index    csub=c{1}(index,:);               % extract the subset of coordinates    else                                % else extract the subset of coordinates and index     csub{1}=c{1}(index,:);    csub{2}=c{2}(index,:);  end;end;