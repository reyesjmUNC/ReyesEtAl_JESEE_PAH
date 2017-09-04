function [] = plotPMvPAH()
% this function will plot PM v all the PAHs

% load data
load('matfiles/pah_data.mat');

% loop through each pah
for i = 1:11    
    idx = ~isnan(val(:,i+4));    
    figure; hold on;
    plot(val(idx,4),val(idx,i+4),'bo');
    xlabel(sprintf('PM_{2.5} (\\mug/m^{3})'));
    ylabel(sprintf('%s (ng/m^{3})',valdispname{i}));
    title(sprintf('PM v %s',valdispname{i}));
    
    % save figure 
    set(gcf,'Position',[0 0 800 500]); 
    print(gcf,'-painters','-dpng','-r600',sprintf('figures/PMvPAH_%s.png',valname{i+4}));
end

end