function [] = readData()
% read in data

[val,valname,filetitle]=readGeoEAS('HardDataNormal.txt');
val(val==-9999) = NaN;

% variables
Longitude = val(:,1);
Latitude = val(:,2);
Time = val(:,3) + datenum([2003 12 31]);
PM2p5 = val(:,4);
benz_a_anthracene = val(:,5);
chrysene = val(:,6);
benzo_bj_fluoranthrene = val(:,7);   
benzo_k_fluoranthrene = val(:,8);
benzo_e_pyrene = val(:,9);
benzo_a_pyrene = val(:,10);
indeno_123cd_pyrene = val(:,11);
benzo_ghi_perylene = val(:,12);
dibenzo_ah_anthracene = val(:,13);
Total_PAHs = val(:,14);
Burden = val(:,15);

timevec = datevec(Time);
timeyr = timevec(:,1);
timemo = timevec(:,2);
timeda = timevec(:,3);

valdispname = {'benz(a)anthracene','chrysene','benzo(bj)fluoranthrene', ...  
    'benzo(k)fluoranthrene','benzo(e)pyrene','benzo(a)pyrene', ...
    'indeno(123cd)pyrene', 'benzo(ghi)perylene','dibenzo(ah)anthracene', ... 
    'Total PAHs','Burden'};

% project longtiude/latitude
cd ../09_mfiles_projections
projxy = ell2lambertcc([Longitude,Latitude],'whiproj2001');
cd ../20_mfiles_PAH
ProjectX = projxy(:,1);
ProjectY = projxy(:,2);
 
% save variables
save('matfiles/pah_data.mat','Longitude','Latitude','Time','PM2p5', ...
    'benz_a_anthracene','chrysene','benzo_bj_fluoranthrene','benzo_k_fluoranthrene', ...
    'benzo_e_pyrene','benzo_a_pyrene','indeno_123cd_pyrene','benzo_ghi_perylene', ...
    'dibenzo_ah_anthracene','Total_PAHs','Burden','Time','timevec','timeyr','timemo','timeda', ...
    'valname','valdispname','val','ProjectX','ProjectY'); 

end