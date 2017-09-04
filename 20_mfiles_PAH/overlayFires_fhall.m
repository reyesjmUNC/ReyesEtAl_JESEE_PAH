function [] = overlayFires_fhall(dataset,spatrad)
% this function will import the fire textfiles and see which estimates lie
% inside or outside the proximity of the fire locations
% Note:
%    hmsfire point proximity, 2/84 hard data in plume
%    hmsfire point proximity, 76/8057 soft data in plume
%    hmsfire point proximity, 298/11726 ck in plume

if nargin < 1, dataset = 1; end
if nargin < 2, spatrad = 50000; end % all data with 50 km from fire locs
% dataset = 1 - overlay hard data 
% dataset = 2 - overlay soft data 
% dataset = 3 - overlay ck

% load data
if dataset == 1
    load('matfiles/pah_data.mat');
    idx = ~isnan(val(:,5));
    cPM = [ProjectX(idx) ProjectY(idx) Time(idx)];
elseif dataset == 2
    load('matfiles/pah_data.mat');
    cPM = [ProjectX ProjectY Time];
else
    load('matfiles/BMEest_krig_PM2p5.mat');
    cPM = ck;
end

% read in fire data files
FID = fopen('PAH_fire/fh_all/firehistory_2005_all_agencies_sub.csv');
datacell = textscan(FID,'%d%s%s%s%f%f%f','HeaderLines',1,'delimiter',',');
fclose(FID);

% get date ranges from start/control date
startnum = cell2mat( cellfun(@datenum,datacell{2},'UniformOutput',false) );
controlnum = cell2mat( cellfun(@datenum,datacell{3},'UniformOutput',false) );

% project longtiude/latitude
longitude = datacell{6};
latitude = datacell{5};
cd ../09_mfiles_projections
projxy = ell2lambertcc([longitude,latitude],'whiproj2001');
cd ../20_mfiles_PAH

% removing all fires less than one acre
idx = datacell{7} >= 1;
startnum = startnum(idx,:);
controlnum = controlnum(idx,:);
projxy = projxy(idx,:);

% reordering fire data
firetime = []; ProjXlong = []; ProjYlong = [];
for i = 1:length(startnum)
    daterange = startnum(i):controlnum(i);
    len = length(daterange);
    firetime = [ firetime ; daterange' ];
    ProjXlong = [ ProjXlong ; repmat(projxy(i,1),len,1) ];
    ProjYlong = [ ProjYlong ; repmat(projxy(i,2),len,1) ];
end
fhall_data = [ProjXlong ProjYlong firetime];

% getting all the days with hard data with hms smoke data date range
unidasmoke = unique(firetime);
idxtime = ismember(cPM(:,3),unidasmoke);
date12 = unique(cPM(idxtime,3));

% loop through each day with both PAH hard data and fh fire data 
inplume = zeros(size(cPM,1),1);
for i = 1:length(date12)

    % subset fire data
    idxfire = firetime == date12(i);
    cFiresub = [ProjXlong(idxfire) ProjYlong(idxfire) firetime(idxfire)];

    % subset relevant PAH locations for day of interest
    cPMsub = cPM(cPM(:,3)==date12(i),1:2);
    
    % all the distances between the fires and hard data
    X = cPMsub; Y = cFiresub(:,1:2);
    DMS = sqrt(bsxfun(@plus,dot(X,X,2),dot(Y,Y,2)')-2*(X*Y'));
    
    % across each row, are there any points within 100 km of a fire?
    firesperpoint = sum(DMS<=spatrad,2);
    inchsub = firesperpoint > 0;
    
    % space/time locations of all points that are in the smoke polygon
    cPMsubtemp = [cPMsub(:,1:2) repmat(date12(i),size(cPMsub,1),1)];
    cPMsubsub = cPMsubtemp(inchsub>0,:);
    
    % putting inplume information back in the final 'inplume' variable
    [aidx bidx] = ismember(cPM,cPMsubsub,'rows');
    inplume(aidx) = 1;
        
end

% turn plume into logical
inplume = logical(inplume);

% save results
disp(nansum(inplume));
if dataset == 1
    save('matfiles/inplume_fhall_hardlocations.mat','inplume','cPM','fhall_data');
elseif dataset == 2
    save('matfiles/inplume_fhall_softlocations.mat','inplume','cPM','fhall_data');    
else
    save('matfiles/inplume_fhall_cklocations.mat','inplume','cPM','fhall_data');
end

end