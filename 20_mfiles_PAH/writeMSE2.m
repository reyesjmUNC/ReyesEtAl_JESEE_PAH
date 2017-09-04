function [] = writeMSE2()
% this will write the MSE results from the cross validation

% load data
load('matfiles/pah_data.mat');
 
fid = fopen('figures/CrossValidation2.txt','wt');
fprintf(fid,'PAH\tsoft MF\tsoftLR\tsoft RAMP\thard\tcokrig\n');

for i = 1:11 % each PAH
    
    fprintf(fid,'%s\t',valname{i+4});
    
    load(sprintf('matfiles/BMEXVal_%s_ellip_MF_soft.mat',valname{i+4}));
    zkall = cell2mat(zkX); zhall = cell2mat(zhX);
    logMSE = mean( (zkall - zhall).^2 );
    fprintf(fid,'%f\t',logMSE);
    
    load(sprintf('matfiles/BMEXVal_%s_ellip_LR_soft.mat',valname{i+4}));
    zkall = cell2mat(zkX); zhall = cell2mat(zhX);
    logMSE = mean( (zkall - zhall).^2 );
    fprintf(fid,'%f\t',logMSE);

    load(sprintf('matfiles/BMEXVal_%s_CAMP.mat',valname{i+4}));
    zkall = cell2mat(zkX); zhall = cell2mat(zhX);
    logMSE = mean( (zkall - zhall).^2 );
    fprintf(fid,'%f\t',logMSE);

    load(sprintf('matfiles/BMEXVal_%s_ellip_NA_hard.mat',valname{i+4}));
    zkall = cell2mat(zkX); zhall = cell2mat(zhX);
    logMSE = mean( (zkall - zhall).^2 );
    fprintf(fid,'%f\t',logMSE);

    load(sprintf('matfiles/BMEXVal_%s_cokrig.mat',valname{i+4}));
    zkall = cell2mat(zkX); zhall = cell2mat(zhX);
    logMSE = mean( (zkall - zhall).^2 );
    fprintf(fid,'%f\n',logMSE);
    
end

fclose(fid);

end