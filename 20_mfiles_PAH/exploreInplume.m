function [] = exploreInplume(pah,method)
% looking at the difference in PAH concentration inside versus outside the
% plume/fire locations at the ck estimation locations
%
% this preforms a right tailed two sample t-test between means assuming 
% unequal variance, seeing if (pah inplume)  > (pah ~inplume)
%
% for a given pah/method, look at all four fire data sources 
%
% method = 1 - soft MF, method = 2 - soft LR, method = 3 - soft CAMP,
% method = 4 - kriging, method = 5 - cokriging
%
% combinations to look at: 11 pahs, 5 methods, 4 fire data sources, 
% buffer size around fire locations (right now, fixed at 50 km)
%

if nargin < 1, pah = 1; end % pah to estimate
if nargin < 2, method = 1; end % estimation method to load

% defining strings
if method == 1
    methstr1 = 'ellip_MF_soft'; methstr2 = 'soft MF';
elseif method == 2
    methstr1 = 'ellip_LR_soft'; methstr2 = 'soft LR';
elseif method == 3
    methstr1 = 'CAMP'; methstr2 = 'soft CAMP';
elseif method == 4
    methstr1 = 'ellip_NA_hard'; methstr2 = 'kriging';
else
    methstr1 = 'cokrig'; methstr2 = 'cokriging';
end

% load data: look at 'valname', 'valdispname'
load('matfiles/pah_data.mat');
    
% load pah prediction data: look at 'ck', 'zk'
load(sprintf('matfiles/BMEest_%s_%s.mat',valname{pah+4},methstr1));

% loop over each fire/smoke data source available
% i = 1 - hmssmoke, i = 2 - hmsfire, i = 3 - fhall, i = 4 - sit
for i = 1:4
    
    % picking fire/smoke data source    
    if i == 1
        firestr = 'hmssmoke';
    elseif i == 2
        firestr = 'hmsfire';
    elseif i == 3
        firestr = 'fhall';
    else
        firestr = 'sit';
    end
    
    % load fire source: look at 'cPM', 'inplume'
    load(sprintf('matfiles/inplume_%s_cklocations.mat',firestr));

    % subsetting prediction to relevent days
    if method == 5 % adding cokriging exception
        [aidx bidx] = ismember(round(ck{1}),round(cPM),'rows');
    else
        [aidx bidx] = ismember(round(ck),round(cPM),'rows');
    end
    bidx(bidx==0) = [];
    zkcompare = zk(bidx);

    % perform t-test
    % h = 1 means the null hypothesis (i.e. the means are equal) is rejected
    if sum(inplume) > 0 
        [zkin lambin] = boxcox(exp(zkcompare(inplume)));
        [zkout lambout] = boxcox(exp(zkcompare(~inplume)));
        [h,p,ci,stats] = ttest2(exp(zkcompare(inplume)),exp(zkcompare(~inplume)),'Vartype','unequal');
        [hin,pin,ksstatin,cvin] = kstest(zkin);
        [hin,pin,ksstatin,cvin] = kstest(zkcompare(inplume));
        [hout,pout,ksstatout,cvout] = kstest(zkcompare(~inplume));
    else
        zkin = NaN; lambin = NaN; zkout = NaN ; lambout = NaN;
        h = NaN; p = NaN; ci = NaN; stats = NaN;
        hin = NaN; pin = NaN; ksstatin = NaN; cin = NaN;
        hout = NaN; pout = NaN; ksstatout = NaN; cout = NaN;
    end

    % save results
    save(sprintf('matfiles/inplume_%s_%s_%s_%s.mat',valname{pah+4}, ...
        methstr1,firestr,'ck'),'cPM','zkcompare','inplume', ...
        'methstr1','methstr2','firestr','h','p','ci','stats', ...
        'hin','pin','ksstatin','cvin','hout','pout','ksstatout','cvout', ...
        'zkin','lambin','zkout','lambout');

    % plot results
    [counts1, values1] = hist(exp(zkcompare(inplume)),100);
    [counts2, values2] = hist(exp(zkcompare(~inplume)),100);
    figure; hold on;
    plot(values1, counts1, 'r-');
    plot(values2, counts2, 'b-');
    legend('inplume','outplume');
    title(sprintf('pah:%s, method:%s, fire source:%s\npredictions: %s, buffsize: %d km, difference: %d, Not Normal in: %d, Not Normal out: %d', ...
        valdispname{pah},methstr2,firestr,'ck',50,h,hin,hout));

    % save figure 
    set(gcf,'Position',[0 0 800 500]); 
    print(gcf,'-painters','-dpng','-r600',sprintf('figures/Histogram_%s_%s_%s_%s.png', ...
        valname{pah+4},methstr1,firestr,'ck'));

end

end