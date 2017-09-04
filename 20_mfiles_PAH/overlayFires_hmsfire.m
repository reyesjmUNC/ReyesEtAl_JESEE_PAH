function [] = overlayFires_hmsfire(dataset,spatrad)
% this function will import the fire textfiles and see which estimates lie
% inside or outside the proximity of the fire locations
% Note:
%    hmsfire point proximity, 30/84 hard data in plume
%    hmsfire point proximity, 1100/8057 soft data in plume
%    hmsfire point proximity, 2261/11726 ck in plume

if nargin < 1, dataset = 3; end
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

% getting all unique days with fire data
allFiles = dir( 'PAH_fire/hmsfire' );
smokefiles = {allFiles(~[allFiles.isdir]).name};
smokeyr = cell2mat( cellfun(@(x) str2double(x(4:7)),smokefiles,'UniformOutput',false) )';
smokemo = cell2mat( cellfun(@(x) str2double(x(8:9)),smokefiles,'UniformOutput',false) )';
smokeda = cell2mat( cellfun(@(x) str2double(x(10:11)),smokefiles,'UniformOutput',false) )';
unidasmoke = unique(datenum([smokeyr smokemo smokeda]));

% getting all the days with hard data with hms fire data date range
idxtime = ismember(cPM(:,3),unidasmoke);
date12 = unique(cPM(idxtime,3));

% loop through each day with both PAH hard data and hms fire data
inplume = zeros(size(cPM,1),1);
for i = 1:length(date12)

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
    cPMsub = cPM(cPM(:,3)==date12(i),1:2);
    
    % all the distances between the fires and hard data
    X = cPMsub; Y = projxy;
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
    save('matfiles/inplume_hmsfire_hardlocations.mat','inplume','cPM');
elseif dataset == 2
    save('matfiles/inplume_hmsfire_softlocations.mat','inplume','cPM');    
else
    save('matfiles/inplume_hmsfire_cklocations.mat','inplume','cPM');
end

end