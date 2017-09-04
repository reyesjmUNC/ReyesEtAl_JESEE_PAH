function [] = writeMSE3()
% this will write the MSE results from the cross validation

% load data
load('matfiles/pah_data.mat');
 
fid = fopen('figures/CrossValidation3.txt','wt');
fprintf(fid,'PAH\tsoft MF\tsoftLR\tsoft RAMP\thard\tcokrig\n');

for i = 1:11 % each PAH

    load(sprintf('matfiles/BMEXVal_%s_ellip_MF_soft.mat',valname{i+4}));
    zkall = cell2mat(zkX); zhall = cell2mat(zhX);
    a = mean( (zkall - zhall).^2 );
        
    load(sprintf('matfiles/BMEXVal_%s_ellip_LR_soft.mat',valname{i+4}));
    zkall = cell2mat(zkX); zhall = cell2mat(zhX);
    b = mean( (zkall - zhall).^2 );

    load(sprintf('matfiles/BMEXVal_%s_CAMP.mat',valname{i+4}));
    zkall = cell2mat(zkX); zhall = cell2mat(zhX);
    c = mean( (zkall - zhall).^2 );

    load(sprintf('matfiles/BMEXVal_%s_ellip_NA_hard.mat',valname{i+4}));
    zkall = cell2mat(zkX); zhall = cell2mat(zhX);
    d = mean( (zkall - zhall).^2 );

    load(sprintf('matfiles/BMEXVal_%s_cokrig.mat',valname{i+4}));
    zkall = cell2mat(zkX); zhall = cell2mat(zhX);
    e = mean( (zkall - zhall).^2 );
    
    [temp ranked] = sort([a b c d e]);
    fprintf(fid,'%s\t%d\t%d\t%d\t%d\t%d\n', ...
        valname{i+4},ranked(1),ranked(2),ranked(3),ranked(4),ranked(5));
    
end

fclose(fid);

end