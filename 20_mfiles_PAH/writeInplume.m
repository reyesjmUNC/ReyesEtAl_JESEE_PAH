function [] = writeInplume()
% this function will collect all the results from the inplume/outofplume
% analysis and put it in a text file

% load data
load('matfiles/pah_data.mat');

% define fixed variables
ckstr = 'ck';
buffstr = '50km';

% initialize text file
fid = fopen('figures/inplume.txt','wt');
fprintf(fid,'PAH\tmethod\tfire\tpredictions\tbuffer\tsignificant\tlb\tub\n');

% loop through all combinations
for i = 1:11 % each pah
    for j = 1:5 % each method
        
        % defining strings
        if j == 1
            methstr1 = 'ellip_MF_soft'; methstr2 = 'soft MF';
        elseif j == 2
            methstr1 = 'ellip_LR_soft'; methstr2 = 'soft LR';
        elseif j == 3
            methstr1 = 'CAMP'; methstr2 = 'soft CAMP';
        elseif j == 4
            methstr1 = 'ellip_NA_hard'; methstr2 = 'kriging';
        else
            methstr1 = 'cokrig'; methstr2 = 'cokriging';
        end

        for k = 1:4 % each fire data set
            
            % picking fire/smoke data source    
            if k == 1
                firestr = 'hmssmoke';
            elseif k == 2
                firestr = 'hmsfire';
            elseif k == 3
                firestr = 'fhall';
            else
                firestr = 'sit';
            end
            
            % load file, add to text file
            load(sprintf('matfiles/inplume_%s_%s_%s_%s.mat',valname{i+4},methstr1,firestr,'ck'));
            fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%d\t%f\t%f\n', ...
                valname{i+4},methstr2,firestr,ckstr,buffstr,h,ci(1),ci(2));
            
        end
    end
end

fclose(fid);

end