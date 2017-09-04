function [] = writeMSE()
% this will write the MSE results from the cross validation

% load data
load('matfiles/pah_data.mat');
 
fid = fopen('figures/CrossValidation.txt','wt');
fprintf(fid, 'PAH\thard/soft\tLR/MF\tlogMSE\tMSE\n');

for i = 1:11 % each PAH
    
    %%% soft data
    
    for j = 0:1 % LR/MF
        
        if j == 1, LRstr = 'LR'; else LRstr = 'MF'; end
        hardstr = 'soft';
        load(sprintf('matfiles/BMEXVal_%s_ellip_%s_%s.mat',valname{i+4},LRstr,hardstr));
        
        zhall = cell2mat(zhX);
        zkall = cell2mat(zkX);
        logMSE = mean( (zkall - zhall).^2 );
        MSE = mean( (exp(zkall) - exp(zhall)).^2 );
        
        fprintf(fid,'%s\t%s\t%s\t%f\t%f\n',valname{i+4},'soft',LRstr,logMSE,MSE);
        
    end
    
    %%% CAMP
    
    load(sprintf('matfiles/BMEXVal_%s_CAMP.mat',valname{i+4}));
        
    zhall = cell2mat(zhX);
    zkall = cell2mat(zkX);
    logMSE = mean( (zkall - zhall).^2 );
    MSE = mean( (exp(zkall) - exp(zhall)).^2 );

    fprintf(fid,'%s\t%s\t%s\t%f\t%f\n',valname{i+4},'soft','CAMP',logMSE,MSE);
    
    %%% hard data only 
    
    LRstr = 'NA';
    hardstr = 'hard';
    
    load(sprintf('matfiles/BMEXVal_%s_ellip_%s_%s.mat',valname{i+4},LRstr,hardstr));
        
    zhall = cell2mat(zhX);
    zkall = cell2mat(zkX);
    logMSE = mean( (zkall - zhall).^2 );
    MSE = mean( (exp(zkall) - exp(zhall)).^2 );

    fprintf(fid,'%s\t%s\t%s\t%f\t%f\n',valname{i+4},'hard',LRstr,logMSE,MSE);
    
    %%% cokriging
    
    LRstr = 'NA';
    hardstr = 'cokrig';
    
    load(sprintf('matfiles/BMEXVal_%s_cokrig.mat',valname{i+4}));
    
    zhall = cell2mat(zhX);
    zkall = cell2mat(zkX);
    logMSE = mean( (zkall - zhall).^2 );
    MSE = mean( (exp(zkall) - exp(zhall)).^2 );

    fprintf(fid,'%s\t%s\t%s\t%f\t%f\n',valname{i+4},'cokrig',LRstr,logMSE,MSE);
    
end

fclose(fid);

end