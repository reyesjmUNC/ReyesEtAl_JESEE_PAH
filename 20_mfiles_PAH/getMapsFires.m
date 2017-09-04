function [] = getMapsFires(pah,LR,soft,showhard,showsoft)
% GETMAPS visualizes results of all pah/soft data/dimension
% combinations
%
% inputs: pah        - This is the name of the PAH being analyzed. This file
%                      should be in single quotes and have the same name as
%                      listed in the heading of PAHData.txt, ie 'chrysene'.
%                      Default: 'benz_a_anthracene'
%                      Type: string
%         LR         - This variable determines which method is used to
%                      find the soft data. If LR = 1, a linear regression 
%                      method will be used to determine the soft data. If 
%                      LR = 0, a mass fraction apprach will be used to 
%                      determine the soft data.
%         soft       - This variable determines hard or soft data used. 
%                      soft = 1 means soft and soft = 0 means hard.
%
% outputs:  
% There are no outputs.
%

if nargin < 1, pah = 1; end
if nargin < 2, LR = 0; end
if nargin < 3, soft = 1; end
if nargin < 4, showhard = 1; end % show hard on maps
if nargin < 5, showsoft = 1; end % show soft on maps
   
% load fire information
load('matfiles/inplume_fhall_cklocations.mat');

% load data
load('matfiles/pah_data.mat');

% loading estimates
if soft == 0, 
    hardstr = 'hard'; 
    LRstr = 'NA';
else
    hardstr = 'soft';
    if LR == 1, LRstr = 'LR'; else LRstr = 'MF'; end
end
load(sprintf('matfiles/BMEest_%s_ellip_%s_%s.mat',valname{pah+4},LRstr,hardstr));

% picking days to show
chdate = datevec(ch(:,3));
unimo = unique(chdate(:,2));
for i = 1:length(unimo)
    idx = unimo(i) == chdate(:,2);
    dayShow(i) = mode(ch(idx,3));    
end

% find all the days in 2005 that have observed data and active fires
unidate1 = unique(ch(:,3));
unidate2 = unique(fhall_data(:,3));
[aidx,bidx] = ismember(unidate1,unidate2);
dayShow = unidate1(aidx);

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
    
    % overlaying soft data
    if soft == 1 & showsoft == 1       
        Property={'Marker','MarkerSize','MarkerEdgeColor'};
        Value ={'o',5,[0 0 0]};       
        idxH = dayShow(i) == cs(:,3);
        colorplot([cs(idxH,1) cs(idxH,2)],exp(zs(idxH)),'redpink',Property,Value,cax);
    end
    
    % plot fire locations
    idxFire = dayShow(i) == fhall_data(:,3);
    plot(fhall_data(idxFire,1),fhall_data(idxFire,2),'kx','MarkerSize',12,'LineWidth',4);

    % plot radial buffer
    th = 0:pi/50:2*pi;
    r = 100000; % 100 km buffer
    fhall_datasub = fhall_data(idxFire,:);
    for j = 1:size(fhall_datasub,1) 
        x = fhall_datasub(j,1);
        y = fhall_datasub(j,2);
        xunit = r * cos(th) + x;
        yunit = r * sin(th) + y;
        plot(xunit,yunit,'k-','LineWidth',3);
    end
    
    set(gca,'XTickLabel','');
    set(gca,'YTickLabel','');
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    text(1800000,-450000,sprintf('ng/m^{3}'),'FontSize',13,'FontWeight','bold');
    if soft == 1
        title(sprintf('Mean %s ellip %s %s on %s with fire',valdispname{pah}, ...
            LRstr,hardstr,datestr(dayShow(i),2)));
    else
        title(sprintf('Mean %s kriging on %s with fire',valdispname{pah},datestr(dayShow(i),2)));
    end
    
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
    if soft == 1
        print(gcf,'-painters','-dpng','-r600', ...
            sprintf('figures/Fire_BMEest_%s_ellip_%s_%s_%s.png', ...
            valname{pah+4},LRstr,hardstr,datestr(dayShow(i))));
    else
        print(gcf,'-painters','-dpng','-r600', ...
            sprintf('figures/Fire_BMEest_%s_%s_%s_%s.png', ...
            valname{pah+4},LRstr,hardstr,datestr(dayShow(i))));
    end

end
     
end