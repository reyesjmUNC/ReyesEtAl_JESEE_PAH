function [] = plot_buffsize()
% plot graph showing the exploratory analysis of buffer size

% load data
load('matfiles/pah_data.mat');

% loop through methods
for i = 1:5
    
    % defining strings
    if i == 1
        methstr1 = 'ellip_MF_soft'; methstr2 = 'soft MF';
    elseif i == 2
        methstr1 = 'ellip_LR_soft'; methstr2 = 'soft LR';
    elseif i == 3
        methstr1 = 'CAMP'; methstr2 = 'soft CAMP';
    elseif i == 4
        methstr1 = 'ellip_NA_hard'; methstr2 = 'kriging';
    else
        methstr1 = 'cokrig'; methstr2 = 'cokriging';
    end
    
    load(sprintf('matfiles/inplume_buffsize_pahcount_%s.mat',methstr1));
    
    plot_hmsfire{i} = pah_hmsfire;
    plot_fhall{i} = pah_fhall;
    
end

% plot figure 1
figure; hold on;
plotstr = {'bo-';'rs-';'co-';'k*-';'rx-'};
for i = 1:5   
    plot(spatrad,plot_hmsfire{i},plotstr{i});    
end
title('buff rad vs sig pahs for hmsfire');
xlabel('buff rad (km)');
ylabel('num of pahs');
legend('softMF','softLR','softCAMP','krig','cokrig');

% save figure 1
set(gcf,'Position',[0 0 800 500]); 
set(gca,'XTickLabel',get(gca,'XTick')/1000);
print(gcf,'-painters','-dpng','-r600','figures/buffrad_hmsfire.png');

% plot figure 2
figure; hold on;
plotstr = {'bo-';'rs-';'co-';'k*-';'rx-'};
for i = 1:5   
    plot(spatrad,plot_fhall{i},plotstr{i});    
end
title('buff rad vs sig pahs for fhall');
xlabel('buff rad (km)');
ylabel('num of pahs');
legend('softMF','softLR','softCAMP','krig','cokrig');

% save figure 2
set(gcf,'Position',[0 0 800 500]); 
set(gca,'XTickLabel',get(gca,'XTick')/1000);
print(gcf,'-painters','-dpng','-r600','figures/buffrad_fhall.png');

end