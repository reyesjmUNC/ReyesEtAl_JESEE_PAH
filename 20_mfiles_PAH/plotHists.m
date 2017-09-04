function [] = plotHists()
% this function will create histograms for each pah

% load data
load('matfiles/pah_data.mat');

% loop through each pah
for i = 1:11
    idx = ~isnan(val(:,i+4));    
    figure; hold on;
    hist(val(idx,i+4),50);
    title(sprintf('histogram of %s',valdispname{i}));
    
    % save figure 
    set(gcf,'Position',[0 0 800 500]); 
    print(gcf,'-painters','-dpng','-r600',sprintf('figures/Histogram_%s.png',valname{i+4}));
end

end