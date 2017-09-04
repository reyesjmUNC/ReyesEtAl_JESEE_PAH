function [] = overlayFires_buffsize()
% this function will explor the ideal buffer size in the hmsfire and fhall
% data sets on the ck  prediction locations

spatrad = 10000:10000:500000; % buffer size in meters

%%% hsmfire data %%%

% load data
load('matfiles/BMEest_krig_PM2p5.mat');
cPM_hmsfire = ck;

% getting all unique days with fire data
allFiles = dir( 'PAH_fire/hmsfire' );
smokefiles = {allFiles(~[allFiles.isdir]).name};
smokeyr = cell2mat( cellfun(@(x) str2double(x(4:7)),smokefiles,'UniformOutput',false) )';
smokemo = cell2mat( cellfun(@(x) str2double(x(8:9)),smokefiles,'UniformOutput',false) )';
smokeda = cell2mat( cellfun(@(x) str2double(x(10:11)),smokefiles,'UniformOutput',false) )';
unidasmoke = unique(datenum([smokeyr smokemo smokeda]));

% getting all the days with hard data with hms fire data date range
idxtime = ismember(cPM_hmsfire(:,3),unidasmoke);
date12 = unique(cPM_hmsfire(idxtime,3));

% loop through each day with both PAH hard data and hms fire data
inplume_hmsfire = zeros(size(cPM_hmsfire,1),length(spatrad));
for i = 1:length(date12)
    disp(i);
    for j = 1:length(spatrad)
        
        % load text file
        datefire = datevec(date12(i));
        FID = fopen(sprintf('PAH_fire/hmsfire/hms%d%0.2d%0.2d.txt', ...
            datefire(1),datefire(2),datefire(3)));
        datacell = textscan(FID,'%f%f%f%s%s','HeaderLines',1,'delimiter',',');
        fclose(FID);
        % project longtiude/latitude
        cd ../09_mfiles_projections
        projxy = ell2lambertcc([datacell{1},datacell{2}],'whiproj2001');
        cd ../20_mfiles_PAH

        % subset relevant PAH locations for day of interest
        cPMsub = cPM_hmsfire(cPM_hmsfire(:,3)==date12(i),1:2);

        % all the distances between the fires and hard data
        X = cPMsub; Y = projxy;
        DMS = sqrt(bsxfun(@plus,dot(X,X,2),dot(Y,Y,2)')-2*(X*Y'));

        % across each row, are there any points within j km of a fire?
        firesperpoint = sum(DMS<=spatrad(j),2);
        inchsub = firesperpoint > 0;

        % space/time locations of all points that are in the smoke polygon
        cPMsubtemp = [cPMsub(:,1:2) repmat(date12(i),size(cPMsub,1),1)];
        cPMsubsub = cPMsubtemp(inchsub>0,:);

        % putting inplume information back in the final 'inplume' variable
        [aidx bidx] = ismember(cPM_hmsfire,cPMsubsub,'rows');
        inplume_hmsfire(aidx,j) = 1;
    
    end  
end

% turn plume into logical
inplume_hmsfire = logical(inplume_hmsfire);

%%% fhall %%%

% load data
load('matfiles/BMEest_krig_PM2p5.mat');
cPM_fhall = ck;

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
    daterange = startnum(i):controlnum(i)+3; 
    len = length(daterange);
    firetime = [ firetime ; daterange' ];
    ProjXlong = [ ProjXlong ; repmat(projxy(i,1),len,1) ];
    ProjYlong = [ ProjYlong ; repmat(projxy(i,2),len,1) ];
end

% getting all the days with hard data with hms smoke data date range
unidasmoke = unique(firetime);
idxtime = ismember(cPM_fhall(:,3),unidasmoke);
date12 = unique(cPM_fhall(idxtime,3));

% loop through each day with both PAH hard data and fh fire data 
inplume_fhall = zeros(size(cPM_fhall,1),length(spatrad));
for i = 1:length(date12)
    disp(i);
    for j = 1:length(spatrad)

        % subset fire data
        idxfire = firetime == date12(i);
        cFiresub = [ProjXlong(idxfire) ProjYlong(idxfire) firetime(idxfire)];

        % subset relevant PAH locations for day of interest
        cPMsub = cPM_fhall(cPM_fhall(:,3)==date12(i),1:2);

        % all the distances between the fires and hard data
        X = cPMsub; Y = cFiresub(:,1:2);
        DMS = sqrt(bsxfun(@plus,dot(X,X,2),dot(Y,Y,2)')-2*(X*Y'));

        % across each row, are there any points within 100 km of a fire?
        firesperpoint = sum(DMS<=spatrad(j),2);
        inchsub = firesperpoint > 0;

        % space/time locations of all points that are in the smoke polygon
        cPMsubtemp = [cPMsub(:,1:2) repmat(date12(i),size(cPMsub,1),1)];
        cPMsubsub = cPMsubtemp(inchsub>0,:);

        % putting inplume information back in the final 'inplume' variable
        [aidx bidx] = ismember(cPM_fhall,cPMsubsub,'rows');
        inplume_fhall(aidx,j) = 1;
    
    end
end

% turn plume into logical
inplume_fhall = logical(inplume_fhall);

% save results
save('matfiles/inplume_buffsize_plusthreedays.mat','cPM_hmsfire','inplume_hmsfire', ...
    'cPM_fhall','inplume_fhall','spatrad'); % added February 28, 2017

end