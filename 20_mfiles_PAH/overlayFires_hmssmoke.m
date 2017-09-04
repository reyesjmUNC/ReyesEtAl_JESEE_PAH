function [] = overlayFires_hmssmoke(dataset)
% this function will import the fire shapefiles and see which estimates lie
% inside or outside the smoke plumes
% Note:
% 1) data in the hmshysplit folder is not used because it's not clear how 
%    to split up the polygons
% 2) hmssmoke polygon overlay, 0/84 hard data in plume
%    hmssmoke polygon overlay, 27/8057 soft data in plume
%    hmssmoke polygon overlay, 94/11726 ck in plume

if nargin < 1, dataset = 3; end
% dataset = 1 - overlay hard data 
% dataset = 2 - overlay soft data 
% dataset = 3 - overlay ck

% load data
if dataset == 1
    load('matfiles/pah_data.mat');
    idx = ~isnan(val(:,5));
    cPM = [val(idx,1:2) datenum([timeyr(idx) timemo(idx) timeda(idx)])]; % estimation locations
elseif dataset == 2
    load('matfiles/pah_data.mat');
    cPM = [val(:,1:2) datenum([timeyr timemo timeda])]; % estimation locations
else
    load('matfiles/BMEest_krig_PM2p5.mat');
    cd ../09_mfiles_projections
    lonlat = lambertcc2ell(ck(:,1:2),'whiproj2001');
    cd ../20_mfiles_PAH
    cPM = [lonlat ck(:,3)];
end

% getting all unique days with smoke data
allFiles = dir( 'PAH_fire/hmssmoke' );
smokefiles = {allFiles(~[allFiles.isdir]).name};
smokeyr = cell2mat( cellfun(@(x) str2double(x(10:13)),smokefiles,'UniformOutput',false) )';
smokemo = cell2mat( cellfun(@(x) str2double(x(14:15)),smokefiles,'UniformOutput',false) )';
smokeda = cell2mat( cellfun(@(x) str2double(x(16:17)),smokefiles,'UniformOutput',false) )';
unidasmoke = unique(datenum([smokeyr smokemo smokeda]));

% getting all the days with hard data with hms smoke data date range
idxtime = ismember(cPM(:,3),unidasmoke);

% subsetting to days that have both smoke data and hard data
cPM = cPM(idxtime,:);

% loop through each day with both PAH hard data and hms smoke data
date12 = unique(cPM(:,3));
inplume = zeros(size(cPM,1),1);
for i = 1:length(date12)

    % load shape file
    datefire = datevec(date12(i));
    hmssmoke = shaperead(sprintf('PAH_fire/hmssmoke/hms_smoke%d%0.2d%0.2d.shp', ...
        datefire(1),datefire(2),datefire(3)));
    shplen = length(hmssmoke);
    
    % subset relevant PAH locations for day of interest
    cPMsub = cPM(cPM(:,3)==date12(i),1:2);
    
    % loop through each polygon in shp to see which points are in the polygon
    inchsub = zeros(size(cPMsub,1),1);
    for j = 1:shplen
        in = inpolygon(cPMsub(:,1),cPMsub(:,2),hmssmoke(j).X,hmssmoke(j).Y);
        inchsub(in) = 1;
    end
    
    % space/time locations of all points that are in the smoke polygon
    cPMsubtemp = [cPMsub(:,1:2) repmat(date12(i),size(cPMsub,1),1)];
    cPMsubsub = cPMsubtemp(inchsub>0,:);
    
    % putting inplume information back in the final 'inplume' variable
    [aidx bidx] = ismember(cPM,cPMsubsub,'rows');
    inplume(aidx) = 1;
        
end

% project cPM
cd ../09_mfiles_projections
projxy = ell2lambertcc(cPM(:,1:2),'whiproj2001');
cd ../20_mfiles_PAH
cPM = [projxy cPM(:,3)];

% turn plume into logical
inplume = logical(inplume);

% save results
disp(nansum(inplume));
if dataset == 1
    save('matfiles/inplume_hmssmoke_hardlocations.mat','inplume','cPM');
elseif dataset == 2
    save('matfiles/inplume_hmssmoke_softlocations.mat','inplume','cPM');    
else
    save('matfiles/inplume_hmssmoke_cklocations.mat','inplume','cPM');
end

end