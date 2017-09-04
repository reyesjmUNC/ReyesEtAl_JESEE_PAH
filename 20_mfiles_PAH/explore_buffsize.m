function [] = explore_buffsize(method)
% this function will explor the ideal buffer size in the hmsfire and fhall
% data sets on the ck  prediction locations
%
% for each data set/buffer size/method, I will could the number of PAH
% where there is a statistically significant difference in inplume/out of
% plume
%
% method = 1 - soft MF, method = 2 - soft LR, method = 3 - soft CAMP,
% method = 4 - kriging, method = 5 - cokriging

if nargin < 1, method = 1; end % estimation method to load

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

% load results: look at 'cPM_hmsfire', 'inplume_hmsfire', 
% 'cPM_fhall', 'inplume_fhall', 'spatrad'
load('matfiles/inplume_buffsize_plusthreedays.mat'); % added February 28, 2017

% loop over each pah
pah_hmsfire = zeros(length(spatrad),1);
pah_fhall = zeros(length(spatrad),1);
h_hmsfire_all = NaN*ones(11,length(spatrad));
ci_hmsfire_all = cell(11,length(spatrad));
h_fhall_all = NaN*ones(11,length(spatrad));
ci_fhall_all = cell(11,length(spatrad));
for i = 1:11
    disp(i);
    for j = 1:length(spatrad)
        
        %%% hmsfire %%%

        % load pah prediction data: look at 'ck', 'zk'
        load(sprintf('matfiles/BMEest_%s_%s.mat',valname{i+4},methstr1));

        % subsetting prediction to relevent days
        if method == 5 % adding cokriging exception
            [aidx bidx] = ismember(round(ck{1}),round(cPM_hmsfire),'rows');
        else
            [aidx bidx] = ismember(round(ck),round(cPM_hmsfire),'rows');
        end
        bidx(bidx==0) = [];
        zkcompare = zk(bidx);

        % perform t-test
        % h = 1 means the null hypothesis (i.e. the means are equal) is rejected
        if sum(inplume_hmsfire) > 0 
            [h,p,ci,stats] = ttest2(exp(zkcompare(inplume_hmsfire(:,j))),exp(zkcompare(~inplume_hmsfire(:,j))), ...
                'Vartype','unequal');
            h_hmsfire_all(i,j) = h;
            ci_hmsfire_all{i,j} = ci;
            if h == 1, pah_hmsfire(j) = pah_hmsfire(j) + 1; end
        else
            h = NaN; p = NaN; ci = NaN; stats = NaN;
        end

        %%% fhall %%%

        % load pah prediction data: look at 'ck', 'zk'
        load(sprintf('matfiles/BMEest_%s_%s.mat',valname{i+4},methstr1));

        % subsetting prediction to relevent days
        if method == 5 % adding cokriging exception
            [aidx bidx] = ismember(round(ck{1}),round(cPM_fhall),'rows');
        else
            [aidx bidx] = ismember(round(ck),round(cPM_fhall),'rows');
        end
        bidx(bidx==0) = [];
        zkcompare = zk(bidx);

        % perform t-test
        % h = 1 means the null hypothesis (i.e. the means are equal) is rejected
        if sum(inplume_fhall) > 0 
            [h,p,ci,stats] = ttest2(exp(zkcompare(inplume_fhall(:,j))),exp(zkcompare(~inplume_fhall(:,j))), ...
                'Vartype','unequal');
            h_fhall_all(i,j) = h;
            ci_fhall_all{i,j} = ci;
            if h == 1, pah_fhall(j) = pah_fhall(j) + 1; end
        else
            h = NaN; p = NaN; ci = NaN; stats = NaN;
        end
    
    end
end

% save results
save(sprintf('matfiles/inplume_buffsize_pahcount_%s_plusthreedays.mat',methstr1), ...
    'spatrad','pah_hmsfire','pah_fhall','h_hmsfire_all','ci_hmsfire_all', ...
    'h_fhall_all','ci_fhall_all'); % added February 28, 2017

end