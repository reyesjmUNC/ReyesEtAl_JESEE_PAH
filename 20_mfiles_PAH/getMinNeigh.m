function [] = getMinNeigh()
% GETMINNEIGH takes the results of analysis.m to find the best
% neighborhood for each PAH, LR, cylinder, and goodNeigh combination

% load data
load('matfiles/pah_data.mat');

% create text file
fid = fopen('figures/optimumSoftDataDevelopment.txt','wt');
fprintf(fid, 'PAH\tCylin\tLR\tminR\tminC\tminMSE\n');

for i = 1:11 % loop through PAH
    
    pah = valname{i+4};
    
    for j = 1:2 % loop through cylin/ellip
        for k = 1:2 % loop through LR/MF

            % defining strings
            if k == 1, LRStr = 'LR'; else LRStr = 'MF'; end
            if j == 1, cyStr = 'cylin'; else cyStr = 'ellip'; end

            % loading data
            load(sprintf('matfiles/exhaustive_%s_%s.mat',pah,cyStr));

            % finding minimum MSE
            if k == 1
                [row col] = find(min(MSELR(:))==MSELR);
                row = row(1); col = col(1); % in case the minimum is not unique
                minMSE = min(MSELR(:));
            else
                [row col] = find(min(MSEMF(:))==MSEMF);
                row = row(1); col = col(1); % in case the minimum is not unique
                minMSE = min(MSEMF(:));
            end

            if j == 1
                minR = spatrad(row);
                minC = temprad(col);
            else
                minR = nmax(row);
                minC = stmetric(col);
            end

            save(sprintf('matfiles/bestMinNeigh_%s_%s_%s.mat',pah,cyStr,LRStr), ...
                'minMSE','minR','minC');
            
             fprintf(fid, '%s\t%s\t%s\t%f\t%f\t%f\n',valname{i+4},cyStr,LRStr,minR,minC,minMSE );  

        end
    end
end

fclose(fid);

end