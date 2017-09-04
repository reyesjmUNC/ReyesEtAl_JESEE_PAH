function [] = getMaps_PM(showhard)
% GETMAPS_PM visualizes results of all PM2.5 data

if nargin < 1, showhard = 1; end % show hard data in plot

% load data
load('matfiles/pah_data.mat');
load(sprintf('matfiles/BMEest_%s_ellip_%s_%s.mat',valname{5},'MF','soft'));

% picking days to show
chdate = datevec(ch(:,3));
unimo = unique(chdate(chdate(:,1)==2005,2));
for i = 1:length(unimo)
    idx = unimo(i) == chdate(:,2);
    dayShow(i) = mode(ch(idx,3));    
end

% loading estimates
load('matfiles/BMEest_krig_PM2p5.mat')

% get state outline
load('../09_mfiles_projections/USAstates5.mat');
cd ../09_mfiles_projections
maskcontour = [ X{39} Y{39} ];
plotax = ell2lambertcc(maskcontour,'whiproj2001');
cd ../20_mfiles_PAH

% loop through each day in a month
for i = 1:length(dayShow)
    
    % picking day
    idx = dayShow(i) == ck(:,3);

    % plotting mean
    ax = [1133721 1940155 -583367 -78441]; % hard coded
    lax = [min(plotax(:,1))-100000 max(plotax(:,1))+100000 min(plotax(:,2))-100000 max(plotax(:,2))+100000];
    colorz = redpink;
    [xg yg Zg] = plotField(ck(idx,:),zkEXP(idx),lax,[plotax(:,1) plotax(:,2)],colorz);
    cax = [prctile(exp(zh),5) prctile(exp(zh),95)];
    caxis(cax);
    colorbar;
    axis(ax);

    % overlaying the counties
    % obtained from: https://connect.ncdot.gov/resources/gis/pages/gis-data-layers.aspx
    % then converted from NC plane to WGS84 via ArcGIS
    nc_counties = shaperead('CountyBoundary_Project/CountyBoundary_Project.shp');
    len = length(nc_counties);
    cd ../09_mfiles_projections
    for j = 1:len        
        counties = ell2lambertcc([nc_counties(j).X',nc_counties(j).Y'],'whiproj2001');        
        plot(counties(:,1),counties(:,2),'k-');
    end
    cd ../20_mfiles_PAH

    % overlaying hard data
    if showhard == 1       
        Property={'Marker','MarkerSize','MarkerEdgeColor'};
        Value ={'s',5,[0 0 0]};       
        idxH = dayShow(i) == ch(:,3);
        colorplot([ch(idxH,1) ch(idxH,2)],exp(zh(idxH)),'redpink',Property,Value,cax);
    end
    
    set(gca,'XTickLabel','');
    set(gca,'YTickLabel','');
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    text(1800000,-450000,sprintf('\\mug/m^{3}'),'FontSize',13,'FontWeight','bold');
    title(sprintf('PM2.5 kriging on %s',datestr(dayShow(i),2)));

    % getting rid of the some the white space surrounding figure
    daspect(gca,[1 1 1]);
    ti = get(gca,'TightInset');
    set(gca,'Position',[ti(1) ti(2)-0.1 1-ti(3)-ti(1) 1-ti(4)-ti(2)+0.2]);
    set(gca,'units','centimeters');
    pos = get(gca,'Position');
    ti = get(gca,'TightInset');

    % save figure 
    set(gcf,'Position',[0 0 800 500]); 
    set(gcf,'PaperUnits','centimeters');
    set(gcf,'PaperSize', [pos(3)+ti(1)+ti(3)+2 pos(4)+ti(2)+ti(4)]+0.2);
    set(gcf,'PaperPositionMode', 'manual');
    set(gcf,'PaperPosition',[0 0 pos(3)+ti(1)+ti(3)+2 pos(4)+ti(2)+ti(4)-2.2]);
    print(gcf,'-painters','-dpng','-r600', ...
        sprintf('figures/BMEest_krig_PM2p5_%s.png',datestr(dayShow(i))));

end
     
end