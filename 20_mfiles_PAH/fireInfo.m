function [] = fireInfo()
% this function will explore 1) which days in 2005 had the most amount of
% fires above 1 acre and 2) which days had the most acres burned (note: I
% only know the total acreage burned during the entirety of the fire and
% NOT which days within that time period burned the most acres)

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

% get size of fire
firesize = datacell{7};

% removing all fires less than one acre
idx = datacell{7} >= 1;
startnum = startnum(idx,:);
controlnum = controlnum(idx,:);
projxy = projxy(idx,:);
firesize = firesize(idx,:);

% reordering fire data
firetime = []; ProjXlong = []; ProjYlong = []; fireSize = [];
acreperday = NaN*ones(length(firesize),1);
daysinrange = cell(length(firesize),1);
for i = 1:length(startnum)
    daterange = startnum(i):controlnum(i);
    len = length(daterange);
    acreperday(i) = firesize(i)/len;
    daysinrange{i} = daterange;
    firetime = [ firetime ; daterange' ];
    ProjXlong = [ ProjXlong ; repmat(projxy(i,1),len,1) ];
    ProjYlong = [ ProjYlong ; repmat(projxy(i,2),len,1) ];
    fireSize = [fireSize; repmat(firesize(i),len,1) ];
end

% intersect all fire days with days that have observed data
load(sprintf('matfiles/BMEest_%s_ellip_%s_%s.mat','benz_a_anthracene','MF','soft'));
unidays = unique(ck(:,3));

idx1 = cellfun(@(x) intersect(unidays,x),daysinrange,'UniformOutput',false);
idx2 = ~cell2mat(cellfun(@isempty,idx1,'UniformOutput',false)); % id of which fires overlap with obs data
acreperdayobs = acreperday(idx2);

% March 5, 2005 is in the 90th percentile of average acres burned/day 
% amoung fires that burned during days with observed PAH data in 2005

daysobsfire = intersect(unidays,firetime);
idx = ismember(firetime,daysobsfire) & fireSize>50;
firetimeobs = firetime(idx); fireSizeobs = fireSize(idx);
unidayfireobs = unique(firetimeobs);

save('firetimeobs.mat','unidayfireobs'); % save these fire/obs days

% find the day that has the most acres burned per day (e.g. acres divided
% the length of time between the start day and control day)

% get count across unique days
a = unique(firetime);
out = [a,histc(firetime(:),a)];
figure; hold on;
plot(out(:,1),out(:,2),'b.');
disp([datevec(out(out(:,2)>=5,1)) out(out(:,2)>=5,2) ])% display days with more than 5 fires
largenumfires = datevec(out(out(:,2)>=5,1));

% get cummulative sum of fire size across unique days
[a,~,c] = unique(firetime);
out = [a, accumarray(c,fireSize)];
figure; hold on;
plot(out(:,1),out(:,2),'b.');
disp([datevec(out(out(:,2)>=300,1)) out(out(:,2)>=300,2) ])% display days with more than 300 acres burned
largeacresfires = datevec(out(out(:,2)>=300,1));

% collect all the day that have a >=5 fires and/or >=300 acres burned and
% pick a day to show as the ratio
largefires = datenum(unique([largenumfires;largeacresfires],'rows'));
largefires = [732376 ; 732406 ; 732418 ; 732577 ; 732630]; % intersected with observed data
save('largefires.mat','largefires')% save these largefires

end