function [] = plot_CI()
% this function will look at the confidence intervals by spatial radius for
% each pah/data source. Each plot will have 5 CI's - one for each
% estimation method.

% load data
load('matfiles/pah_data.mat');

for i = 1:2 % fhall and hmsfire
    
    if i == 1, firestr = 'fhall'; else firestr = 'hmsfire'; end
    
    for j = 1:11 % each pah  
        
        figure; hold on;
        plotstr = {'bo-','bs-','r*-','b.-','c^-'};
        disp([i,j]);
        
        for k = 1:5 % each method
            
            % defining strings
            if k == 1
                methstr1 = 'ellip_MF_soft'; methstr2 = 'soft MF';
            elseif k == 2
                methstr1 = 'ellip_LR_soft'; methstr2 = 'soft LR';
            elseif k == 3
                methstr1 = 'CAMP'; methstr2 = 'soft CAMP';
            elseif k == 4
                methstr1 = 'ellip_NA_hard'; methstr2 = 'kriging';
            else
                methstr1 = 'cokrig'; methstr2 = 'cokriging';
            end
        
            % load results
            load(sprintf('matfiles/inplume_buffsize_pahcount_%s.mat',methstr1));
            
            % plotting
            if i == 1
                idx = h_fhall_all(j,:) == 1; 
                lb = cellfun(@(v) v(1), ci_fhall_all(j,:));
                ub = cellfun(@(v) v(2), ci_fhall_all(j,:));
            else
                idx = h_hmsfire_all(j,:) == 1; 
                lb = cellfun(@(v) v(1), ci_hmsfire_all(j,:));
                ub = cellfun(@(v) v(2), ci_hmsfire_all(j,:));
            end
            
            % keeping legend
            if k == 1
                h1 = plot(spatrad(idx),lb(idx),plotstr{k},spatrad(idx),ub(idx),plotstr{k});
            elseif k == 2
                h2 = plot(spatrad(idx),lb(idx),plotstr{k},spatrad(idx),ub(idx),plotstr{k});
            elseif k == 3
                h3 = plot(spatrad(idx),lb(idx),plotstr{k},spatrad(idx),ub(idx),plotstr{k});
            elseif k == 4
                h4 = plot(spatrad(idx),lb(idx),plotstr{k},spatrad(idx),ub(idx),plotstr{k});
            else
                h5 = plot(spatrad(idx),lb(idx),plotstr{k},spatrad(idx),ub(idx),plotstr{k});
            end
        
        end

        % save figure        
        a = get(gca,'XLim'); 
        plot(a,[0 0],'k--');
        if ~isempty(h5) & isempty(h3)
            legend([h1(1) h2(1) h4(1) h5(1)],{'soft MF','soft LR','krig','cokrig'});
        elseif ~isempty(h5)
            legend([h1(1) h2(1) h3(1) h4(1) h5(1)],{'soft MF','soft LR','soft CAMP','krig'});
        else
            legend([h1(1) h2(1) h3(1) h4(1)],{'soft MF','soft LR','soft CAMP','krig'});
        end
        xlabel('meters');
        ylabel('confidence interval');
        title(sprintf('%s %s',firestr,valdispname{j}));
        set(gcf,'Position',[0 0 800 500]); 
        print(gcf,'-painters','-dpng','-r600', ...
            sprintf('figures/buffrad_CI_%s_%s.png',firestr,valname{j+4}));
        
    end
    
    close all;
    
end

end