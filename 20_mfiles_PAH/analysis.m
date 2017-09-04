function [] = analysis()
% this function will call all the functions in the PAH directory and
% perform the entire analysis needed for this project

% loading BME function
cd ../BMELIB2.0b
startup;
cd ../20_mfiles_PAH

% reading raw data
readData();

%%%%%%%%%%%%%%%%%%%% STAGE ONE: EXPLORATORY ANALYSIS %%%%%%%%%%%%%%%%%%%%%

% plotting histograms of PAHs
plotHists();
close all;

% plot PM v PAH
plotPMvPAH();
close all;

%%%%%%%%%%%%%%%%%%%% STAGE TWO: SOFT DATA DEVELOPMENT %%%%%%%%%%%%%%%%%%%%

% performing an exhaustive search on the soft data development
% this function needs to be run on the cluster by typing:
% 'sh exhaustiveSoftDevelopment.sh'
exhaustiveSoftDevelopment();
close all;

% visualizing the results of exhaustiveSoftDevelopment.m
visualizeSoftDevelopment();
close all;

% picking best parameters for the soft data development
% (cylindrical/elliposidal and linear regression/mass fraction) and
% displying results in a text file
getMinNeigh();

%%%%%%%%% STAGE THREE: ESTIMATING SOFT DATA AND COVARIANCE MODEL %%%%%%%%%

% find covariance models
for i = 1:11
    for j = 0:1
        disp([i j]);
        autoCovmodel(i,j);
    end
end
autoCovModel_PM();
close all;

% estimating soft data - linear regression and mass fraction method
for i = 1:11
    for j = 0:1
        disp([i j]);
        estSoftData(i,j);
    end
end

%%%%%%%%%%%%%%%%%%%%%% STAGE FOUR: CROSS-VALIDATION %%%%%%%%%%%%%%%%%%%%%%

% performing all different cross-validations 
for i = 1:11 % each PAH
    for j = 0:1 % LR/MF
        disp([i j]);
        getValidation(i,j,1);
    end
end

%%%%%%%%%%%%%%%%%%%%% STAGE FIVE: COMPARISON METHODS %%%%%%%%%%%%%%%%%%%%%

% kriging 
for i = 1:11 % each PAH
    disp(i);
    getValidation(i,0,0);
end

% cokriging runs
for i = 1:11
    getCoKrigVar(i);
end
close all;

for i = 1:11
    disp(i);
    getCoKrigValidation(i);
end

%%%%%%%%%%%%%%%%%%%%%%% STAGE SIX: COMPARE METHODS %%%%%%%%%%%%%%%%%%%%%%%

% write results into a file
writeMSE();
writeMSE2();
writeMSE3();

%%%%%%%%%%%%%%%%%%%%%%%%% STAGE SEVEN: MAKE MAPS %%%%%%%%%%%%%%%%%%%%%%%%%

% running kriging analysis to create maps
getBMEest_PM();

% running BME analysis to create maps
for i = 1:11
    for j = 0:1
        disp([i j]);
        getBMEest(i,j,1);
    end
end

% running kriging analysis to create map
for i = 1:11
    disp(i);
    getBMEest(i,0,0);
end

% running cokriging analysis to create maps
for i = 1:11
    disp(i);
    getBMEest_cokrig(i);
end

% create PM maps
getMaps_PM();
close all;

% create PAH maps from soft data
for i = 1:11
    for j = 0:1
        getMaps(i,j,1);
        close all;
    end
end

% create PAH maps from hard data
for i = 1:11
    getMaps(i,0,0);
    close all;
end

% create PAH maps for cokriging data
for i = 1:11
    getMaps_cokrig(i);
    close all;
end

%%%%%%%%%%%%%%%%%%% STAGE EIGHT: EXPLORE PAH AND FIRES %%%%%%%%%%%%%%%%%%%

% tying plume/fire information to PAHs
for i = 1:3
    disp(i);
    overlayFires_hmssmoke(i);
    overlayFires_hmsfire(i);
    overlayFires_fhall(i);
    overlayFires_sit(i);
end

% explore difference inside/outside fires
for i = 1:11 % each pah
    for j = 1:5 % each estimation method
        disp([i,j]);
        exploreInplume(i,j);
        close all;
    end
end

% gather all explore plume methods into one file
writeInplume();
% conclusions: exclude hmssmoke and sit fire sources

% exploring km buffer size
overlayFires_buffsize();
for i = 1:5
    explore_buffsize(i);
end
plot_buffsize();
% conclusions: 100 km for fhall

% explore conf intervals of differences in means of inplume/outofplume
plot_CI();
% conclusions: 100 km for fhall

% explore fire info: days with the most fires, days with the largest fires
fireInfo(); % 11/18/2016 KEEP WORKING ON THIS

% overlay PAH maps with the most compelling fire information
for i = 1:11
    for j = 0:1
        getMapsFires(i,j,1);
        close all;
    end
end

% create PAH maps from hard data
for i = 1:11
    getMapsFires(i,0,0);
    close all;
end

% create PAH maps for cokriging data
for i = 1:11
    getMapsFires_cokrig(i);
    close all;
end

end