function [] = visualizeSoftDevelopment()
%
% VISUALIZESOFTEVELOPMENT looks at MSE at a function of PAH and (LR cylinder or
% MF cylinder or LR egg or MF egg). This will also compare the old and new
% approach. The old approach forces all neighborhoods to have a minimum
% number of neighbors. The new approach does not have this restriction. The
% goal of this function is to determine if the old/new method is better and
% which neighborhood is better.

% load data
load('matfiles/pah_data.mat');

for i = 1:11 % loop through each PAH
    
    pah = valname{i+4};
    
    for j = 1:2 % loop through cylinder/ellip
        for k = 1:2 % loop through LR/MF

            % defining strings
            if j == 1 
                cyStr = 'cylin'; 
                cyLabel = 'cylindrical';
                xLabel = 'days';
                yLabel = 'kilometers';
            else
                cyStr = 'ellip';
                cyLabel = 'ellipsoidal';
                xLabel = 'space/time metric';
                yLabel = 'number of points';
            end
            if k == 1
                LRLabel = 'linear regression';
                LRStr = 'LR';
            else
                LRLabel = 'mass fraction';
                LRStr = 'MF';
            end

            % loading data
            load(sprintf('matfiles/exhaustive_%s_%s.mat',pah,cyStr));

            % plotting data
            figure; hold on;
            colormap(flipud('redpink(16)')); % inverses the color order

            % marking where the minimum is
            if k == 1
                [xloc yloc] = find(MSELR==min(MSELR(:))); xloc=xloc(1); yloc=yloc(1);
                if j == 1
                    markX = temprad(yloc); markY = spatrad(xloc);
                else
                    markX = stmetric(yloc); markY = nmax(xloc);
                end
            else
                [xloc yloc] = find(MSEMF==min(MSEMF(:))); xloc=xloc(1); yloc=yloc(1);
                if j == 1
                    markX = temprad(yloc); markY = spatrad(xloc);
                else
                    markX = stmetric(yloc); markY = nmax(xloc);
                end
            end

            if j == 1
                [x y] = meshgrid(temprad,spatrad);
                if k == 1
                    pcolor(x,y,MSELR);
                else
                    pcolor(x,y,MSEMF);
                end
            else
                [x y] = meshgrid(stmetric,nmax);
                if k == 1
                    pcolor(x,y,MSELR);
                else
                    pcolor(x,y,MSEMF);
                end
            end

            % marking the minimum MSE
            plot(markX,markY,'gx','MarkerSize',10,'LineWidth',2);

            if j == 1
                idx = ~isnan(MSELR); minx = x(idx);
                idx = ~isnan(MSELR); miny = y(idx);
            else
                axis([min(x(:)) max(x(:)) min(y(:)) max(y(:))]);
            end
            
            % get color axis
            caxis([0 2.5]);
            
            shading interp;
            colorbar;
            xlabel(sprintf('%s',xLabel));
            ylabel(sprintf('%s',yLabel));
            if j == 1, set(gca,'YTickLabel',get(gca,'YTick')./1000); end
            title(sprintf('MSE %s for an %s neighborhood \nusing %s approach', ...
                valdispname{i},cyLabel,LRLabel));

            % save figure 
            set(gcf,'Position',[0 0 800 500]); 
            print(gcf,'-painters','-dpng','-r600',sprintf('figures/exhaustive_%s_%s_%s.png',pah,cyStr,LRStr));

        end
    end
end

end