%%
maindir = '/Users/psyc1586_admin/GVB_data/harmonisation_work/00_aa_skullImprovement_Oct2025/00_running_for_paper/run_eval_metrics_IQMscannerOnly/data';
fname   = 'results.mat';
load(fullfile(maindir, fname));

ori_idps  =     {'T1_FIRST_left_hippocampus',...
                'T1_FIRST_right_hippocampus', ...
                'T1_SIENAX_CSF_norm_vol', ...
                'T1_SIENAX_GM_norm_vol', ...
                'T1_SIENAX_WM_norm_vol', ...
                'T1_SIENAX_brain_norm_vol', ...
                'T1_SIENAX_periphGM_norm_vol'};
plot_idps = {'Hipp-L','Hipp-R','CSF','GM', 'WM','Brain', 'PeriphGM'};
idpMap    = containers.Map(ori_idps, plot_idps);

% for biological variablity use this threshold - bonferroni corrected
p_thr = round(0.05/length(ori_idps),3);

methodOrder = {'Raw', ...
               'HistMatch', ...
               'SynthSR', ...
               'LR-Scanner', ...
               'ComBat-Scanner', ...
               'CovBat-Scanner', ...
               'LME-Scanner', ...
               'LME-IQM-Scanner',...
               'LongComBat-Scanner', ...
               'LME-ScannerFSiteR', ...
               'LongComBat-ScannerFSiteR', ...
               'LR-Site', ...
               'ComBat-Site', ...
               'CovBat-Site', ...
               'LME-Site', ...
               'LME-IQM-Site', ...
               'LongComBat-Site', ...
               };
idpOrder = {'Brain', 'GM', 'PeriphGM', 'WM', 'CSF', 'Hipp-L', 'Hipp-R'};

idp_names      = {'Brain','CSF','GM','Hipp-L','Hipp-R','PeriphGM','WM'};
tissue_markers = {'o', '^', 's', 'v', '>', '<', 'p'};  % 7 distinct shapes
tissue_colors  = [
    0.8500 0.3250 0.0980;  
    0.9290 0.6940 0.1250; 
    0.4940 0.1840 0.5560;  
    0.3010 0.7450 0.9330;  
    0.6350 0.0780 0.1840;  
    0.0000 0.4470 0.7410; 
    1.0000 0.0000 1.0000  
];
markerMap = containers.Map(idp_names, tissue_markers);
colorMap = containers.Map(idp_names, num2cell(tissue_colors,2));


%% residual variability plots
res_var         = allResults.residualvariability;
res_var.idpname = cellfun(@(x) idpMap(x), res_var.idpname, 'UniformOutput', false);

% 1) overall site effect
res_var_site         = res_var(ismember(res_var.Effects, 'site'),:);
res_var_site.method  = categorical(res_var_site.method, methodOrder, 'Ordinal', true);
res_var_site.idpname = categorical(res_var_site.idpname, idpOrder, 'Ordinal', true);

h = heatmap(res_var_site, "method", "idpname","ColorVariable", "anova_batches", ...
            "Colormap", parula);
% Remove title and axis labels
h.Title = '';
h.XLabel = '';
h.YLabel = '';
h.ColorbarVisible='off';
% Adjust font size
h.FontSize = 14;
saveas(gcf, 'residualVariability_overallSite.fig')
print('residualVariability_overallSite', '-dpng', '-r900');
close(gcf)

% 2) pairwise site effect
h = heatmap(res_var_site, "method", "idpname","ColorVariable", "n_is_batchSig", ...
            "Colormap", parula);
% Remove title and axis labels
h.Title = '';
h.XLabel = '';
h.YLabel = '';
h.ColorbarVisible='on';
% Adjust font size
h.FontSize = 14;
saveas(gcf, 'residualVariability_pairwiseSite.fig')
print('residualVariability_pairwiseSite', '-dpng', '-r900');
close(gcf)

% 1) overall scanner effect
res_var_scanner         = res_var(ismember(res_var.Effects, 'scanner'),:);
res_var_scanner.method  = categorical(res_var_scanner.method, methodOrder, 'Ordinal', true);
res_var_scanner.idpname = categorical(res_var_scanner.idpname, idpOrder, 'Ordinal', true);

h = heatmap(res_var_scanner, "method", "idpname","ColorVariable", "anova_batches", ...
            "Colormap", parula);
% Remove title and axis labels
h.Title = '';
h.XLabel = '';
h.YLabel = '';
h.ColorbarVisible='off';
% Adjust font size
h.FontSize = 14;
saveas(gcf, 'residualVariability_overallScanner.fig')
print('residualVariability_overallScanner', '-dpng', '-r900');
close(gcf)

%% Between subject variability
betw_var         = allResults.betSubjVariability;
betw_var.idpname = cellfun(@(x) idpMap(x), betw_var.idpname, 'UniformOutput', false);
% 1) Get values for "site" only
betw_var_site        = betw_var(ismember(betw_var.Effects, 'site'),:);
betw_var_icc         = betw_var_site(:,{'method','idpname','ICC'});
betw_var_icc.method  = categorical(betw_var_icc.method, methodOrder, 'Ordinal', true);
betw_var_icc.idpname = categorical(betw_var_icc.idpname, idpOrder, 'Ordinal', true);
% Convert idpname to cell array of chars if necessary
if ~iscell(betw_var_icc.idpname)
    betw_var_icc.idpname = cellstr(betw_var_icc.idpname);
end
figure;
h = boxplot(betw_var_icc.ICC, betw_var_icc.method); % no default outliers
set(h,'Linew', 1.2)
hold on;
jitterAmount = 0.1;
methods      = categories(betw_var_icc.method);
for i = 1:numel(methods)
    idx = betw_var_icc.method == methods{i};
    
    x = i + (rand(sum(idx),1)-0.5)*jitterAmount; % jitter
    y = betw_var_icc.ICC(idx);
    tissues = betw_var_icc.idpname(idx);
    
    for j = 1:numel(y)
        scatter(x(j), y(j), 100, ...
            'Marker', markerMap(tissues{j}), ...
            'MarkerFaceColor', colorMap(tissues{j}), ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.8);
    end
end
hold off;
ylabel('ICC');
set(gca, 'XTickLabel', methodOrder);
set(gca, 'FontSize', 14);
box off;
grid on;
saveas(gcf, 'betweenVariability_ICC.fig')
print('betweenVariability_ICC', '-dpng', '-r900');
close(gcf)

%% Biological variability: Age
age_biovar           = allResults.bioAgeVariability;
age_biovar.idpname   = cellfun(@(x) idpMap(x), age_biovar.idpname, 'UniformOutput', false);
% 1) Get values for "site" only
age_biovar_site      = age_biovar(ismember(age_biovar.Effects, 'site'),:);
age_biovar_site.method  = categorical(age_biovar_site.method, methodOrder, 'Ordinal', true);
age_biovar_site.idpname = categorical(age_biovar_site.idpname, idpOrder, 'Ordinal', true);
age_biovar_site.pval_Age = round(age_biovar_site.pval_Age, 3);

age_table = age_biovar_site;
if ~iscell(age_table.idpname)
    age_table.idpname = cellstr(age_table.idpname);
end

methods = categories(age_table.method);
numIDPs = numel(idp_names);
offset = linspace(-0.3, 0.3, numIDPs);  % fixed offsets per IDP
f = figure;
f.Position(3:4) = [1400, 600];
hold on;

for i = 1:numel(methods)
    idx = age_table.method == methods{i};
    subT = age_table(idx,:);
    
    for j = 1:numel(idp_names)
        tissue = idp_names{j};
        idx_tissue = strcmp(subT.idpname, tissue);
        if any(idx_tissue)
            est = subT.est_Age(idx_tissue);
            ciL = subT.ciL_Age(idx_tissue);
            ciU = subT.ciU_Age(idx_tissue);
            pval = subT.pval_Age(idx_tissue);
            
            % X position offset for tissue
            x = i + offset(j);
            
            % Marker and color
            baseColor = colorMap(tissue);
            mark = markerMap(tissue);
            
            % Determine significance color
            if (pval < p_thr) && (pval < 0.05)
                faceColor = [1 0 0];      % red
            elseif (pval < 0.05) && (pval > p_thr)
                faceColor = [1 0.5 0];
            else
                faceColor = [0.5 0.5 0.5]; % gray
            end
            
            % CI line
            line([x x], [ciL ciU], 'Color', baseColor, 'LineWidth', 2.5);
            
            % Estimate point
            scatter(x, est, 100, 'Marker', 'diamond', ...
                'MarkerFaceColor', faceColor, ...
                'MarkerEdgeColor', 'k', ...
                'LineWidth', 1.2);
        end
    end
end
hold off;
xlim([0.5 numel(methods)+0.5]);
set(gca, 'XTick', 1:numel(methods), 'XTickLabel', methods);
set(gca, 'FontSize', 14, 'Box', 'on');
xlabel('Method');
ylabel('Age assoication: effect size (β), CI, p-vals');
% title('Age associations across methods and IDPs');
yline(0, '--k', 'LineWidth', 1); box off; grid on; 

% --- Custom legends ---
hold on;
% ---- 1. IDP legend (colored lines only) ----
legendHandles_idp = gobjects(numIDPs,1);
for j = 1:numIDPs
    legendHandles_idp(j) = plot(nan, nan, '-', ...
        'Color', colorMap(idp_names{j}), 'LineWidth', 2);
end
% ---- 2. Significance legend (scatter dots) ----
legendHandles_sig(1) = scatter(nan, nan, 100, 'diamond', ...
    'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', 'k', 'DisplayName', sprintf('p < %.3f',round(p_thr,3)));
legendHandles_sig(2) = scatter(nan, nan, 100, 'diamond', ...
    'MarkerFaceColor', [1 0.5 0], 'MarkerEdgeColor', 'k', 'DisplayName', sprintf('%.3f < p < %.2f',round(p_thr,3),0.05));
legendHandles_sig(3) = scatter(nan, nan, 100, 'diamond', ...
    'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', 'k', 'DisplayName', sprintf('p ≥ %.3f',0.05));
% ---- 3. Combine both legends ----
legend([legendHandles_idp; legendHandles_sig'], ...
    [idp_names, {sprintf('p < %.3f',round(p_thr,3)), sprintf('%.3f < p < %.2f',round(p_thr,3),0.05),sprintf('p ≥ %.2f',0.05)}], ...
    'Location', 'eastoutside', 'FontSize', 12);
hold off;
saveas(gcf, 'bioVariability_Age.fig')
print('bioVariability_Age', '-dpng', '-r900');
close(gcf)

%% Biological variability: Timepoint
time_biovar           = allResults.bioTimeVariability;
time_biovar.idpname   = cellfun(@(x) idpMap(x), time_biovar.idpname, 'UniformOutput', false);
% 1) Get values for "site" only
time_biovar_site      = time_biovar(ismember(time_biovar.Effects, 'site'),:);
time_biovar_site.method  = categorical(time_biovar_site.method, methodOrder, 'Ordinal', true);
time_biovar_site.idpname = categorical(time_biovar_site.idpname, idpOrder, 'Ordinal', true);
time_biovar_site.pval_time = round(time_biovar_site.pval_time, 3);

time_table = time_biovar_site;
if ~iscell(time_table.idpname)
    time_table.idpname = cellstr(time_table.idpname);
end

methods = categories(time_table.method);
numIDPs = numel(idp_names);
offset = linspace(-0.3, 0.3, numIDPs);  % fixed offsets per IDP
f = figure;
f.Position(3:4) = [1400, 600];
hold on;

for i = 1:numel(methods)
    idx = time_table.method == methods{i};
    subT = time_table(idx,:);
    
    for j = 1:numel(idp_names)
        tissue = idp_names{j};
        idx_tissue = strcmp(subT.idpname, tissue);
        if any(idx_tissue)
            est = subT.est_time(idx_tissue);
            ciL = subT.ciL_time(idx_tissue);
            ciU = subT.ciU_time(idx_tissue);
            pval = subT.pval_time(idx_tissue);
            
            % X position offset for tissue
            x = i + offset(j);
            
            % Marker and color
            baseColor = colorMap(tissue);
            mark = markerMap(tissue);
            
            % Determine significance color
            if (pval < p_thr) && (pval < 0.05)
                faceColor = [1 0 0];      % red
            elseif (pval < 0.05) && (pval > p_thr)
                faceColor = [1 0.5 0];
            else
                faceColor = [0.5 0.5 0.5]; % gray
            end
            
            % CI line
            line([x x], [ciL ciU], 'Color', baseColor, 'LineWidth', 2.5);
            
            % Estimate point
            scatter(x, est, 100, 'Marker', 'diamond', ...
                'MarkerFaceColor', faceColor, ...
                'MarkerEdgeColor', 'k', ...
                'LineWidth', 1.2);
        end
    end
end
hold off;
xlim([0.5 numel(methods)+0.5]);
set(gca, 'XTick', 1:numel(methods), 'XTickLabel', methods);
set(gca, 'FontSize', 14, 'Box', 'on');
xlabel('Method');
ylabel('Timepoint assoication: effect size (β), CI, p-vals');
% title('Timepoint associations across methods and IDPs');
yline(0, '--k', 'LineWidth', 1); box off; grid on; 

% --- Custom legends ---
hold on;
% ---- 1. IDP legend (colored lines only) ----
legendHandles_idp = gobjects(numIDPs,1);
for j = 1:numIDPs
    legendHandles_idp(j) = plot(nan, nan, '-', ...
        'Color', colorMap(idp_names{j}), 'LineWidth', 2);
end
% ---- 2. Significance legend (scatter dots) ----
legendHandles_sig(1) = scatter(nan, nan, 100, 'diamond', ...
    'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', 'k', 'DisplayName', sprintf('p < %.3f',round(p_thr,3)));
legendHandles_sig(2) = scatter(nan, nan, 100, 'diamond', ...
    'MarkerFaceColor', [1 0.5 0], 'MarkerEdgeColor', 'k', 'DisplayName', sprintf('%.3f < p < %.2f',round(p_thr,3),0.05));
legendHandles_sig(3) = scatter(nan, nan, 100, 'diamond', ...
    'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', 'k', 'DisplayName', sprintf('p ≥ %.3f',0.05));
% ---- 3. Combine both legends ----
legend([legendHandles_idp; legendHandles_sig'], ...
    [idp_names, {sprintf('p < %.3f',round(p_thr,3)), sprintf('%.3f < p < %.2f',round(p_thr,3),0.05),sprintf('p ≥ %.2f',0.05)}], ...
    'Location', 'eastoutside', 'FontSize', 12);
hold off;
saveas(gcf, 'bioVariability_Time.fig')
print('bioVariability_Time', '-dpng', '-r900');
close(gcf)

%% multivariate variablity plots (MD)
mdVals = allResults.multivariateVariability;
% Stack all method columns into a single "Value" column
mdlong = stack(mdVals, 2:width(mdVals), ...        % stack all columns except 'batch'
    'NewDataVariableName','MDval', ...  % name of numeric variable
    'IndexVariableName','method');      % name of method column
mdlong.method = categorical(strrep(string(mdlong.method), '_', '-'));
mdlong.method = categorical(mdlong.method, methodOrder, 'Ordinal', true);
mdlong.MDval  = round(mdlong.MDval, 1);
mdlong.batch  = categorical(strrep(string(mdlong.batch), 'CAM', 'CAM-GE'));
mdlong.batch  = categorical(strrep(string(mdlong.batch), 'EDI', 'EDI-Siemens'));
mdlong.batch  = categorical(strrep(string(mdlong.batch), 'ICL', 'ICL-GE'));
mdlong.batch  = categorical(strrep(string(mdlong.batch), 'KCL', 'KCL-Siemens'));
mdlong.batch  = categorical(strrep(string(mdlong.batch), 'UCL', 'UCL-Siemens'));
mdlong.batch  = categorical(strrep(string(mdlong.batch), 'MAN', 'MAN-GE'));
mdlong.batch  = categorical(strrep(string(mdlong.batch), 'SHF', 'SHF-GE'));
mdlong.batch  = categorical(strrep(string(mdlong.batch), 'NEW', 'NEW-GE'));
mdlong.batch  = categorical(strrep(string(mdlong.batch), 'average_batch', 'Avg-MD'));
siteOrder = {'Avg-MD','EDI-Siemens', ...
              'KCL-Siemens', ...
              'UCL-Siemens', ...
              'CAM-GE',...
              'ICL-GE', ...
              'MAN-GE', ...
              'NEW-GE', ...
              'SHF-GE'};
mdlong.batch = categorical(mdlong.batch, siteOrder, 'Ordinal', true);
f = figure;
f.Position(3:4) = [1000, 600];
h = heatmap(mdlong, "method", "batch","ColorVariable", "MDval", ...
            "Colormap", parula);
% Remove title and axis labels
h.Title = '';
h.XLabel = '';
h.YLabel = '';
h.ColorbarVisible='on';
% Adjust font size
h.FontSize = 14;
saveas(gcf, 'multivariateVariability_MD.fig')
print('multivariateVariability_MD', '-dpng', '-r900');
close(gcf)

%% Within subject variability - subject order
subjorder_struct       = allResults.withinsubjvariability.subjorder;
allfields              = fieldnames(subjorder_struct);
allTables              = cellfun(@(f) addvars(subjorder_struct.(f), repmat({f}, height(subjorder_struct.(f)), 1), ...
                        'NewVariableNames', 'idpname'), allfields, 'UniformOutput', false);
bigTsubjOrder          = vertcat(allTables{:});
bigTsubjOrder.idpname  = cellfun(@(x) idpMap(x), bigTsubjOrder.idpname, 'UniformOutput', false);
bigTsubjOrder.IDP      = strrep(bigTsubjOrder.IDP, "_", "-");
bigTsubjOrder.method   = bigTsubjOrder.IDP;

% 1) Repeat
repeat_subjorder          = bigTsubjOrder(ismember(bigTsubjOrder.groupname, 'Repeat'),:);
repeat_subjorder.method   = categorical(repeat_subjorder.method, methodOrder, 'Ordinal', true);
repeat_subjorder.idpname  = categorical(repeat_subjorder.idpname, idpOrder, 'Ordinal', true);

% Convert idpname to cell array of chars if necessary
if ~iscell(repeat_subjorder.idpname)
    repeat_subjorder.idpname = cellstr(repeat_subjorder.idpname);
end
f = figure;
f.Position(3:4) = [1600,1200];
tiledlayout(1,3)

nexttile;
h = boxplot(repeat_subjorder.SpearmanRho, repeat_subjorder.method); % no default outliers
set(h,'Linew', 1.2)
hold on;
jitterAmount = 0.1;
methods      = categories(repeat_subjorder.method);
for i = 1:numel(methods)
    idx = repeat_subjorder.method == methods{i};
    
    x = i + (rand(sum(idx),1)-0.5)*jitterAmount; % jitter
    y = repeat_subjorder.SpearmanRho(idx);
    tissues = repeat_subjorder.idpname(idx);
    
    for j = 1:numel(y)
        scatter(x(j), y(j), 100, ...
            'Marker', markerMap(tissues{j}), ...
            'MarkerFaceColor', colorMap(tissues{j}), ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.8);
    end
end
hold off;
ylabel('Spearman Correlation');
ylim([-0.3,1.1]);
set(gca, 'XTickLabel', methodOrder);
title('a) Repeat-scanner')
set(gca, 'FontSize', 14);
box off;
grid on;

% 2) Intra
intra_subjorder          = bigTsubjOrder(ismember(bigTsubjOrder.groupname, 'Intra-data'),:);
intra_subjorder.method   = categorical(intra_subjorder.method, methodOrder, 'Ordinal', true);
intra_subjorder.idpname  = categorical(intra_subjorder.idpname, idpOrder, 'Ordinal', true);
% Convert idpname to cell array of chars if necessary
if ~iscell(intra_subjorder.idpname)
    intra_subjorder.idpname = cellstr(intra_subjorder.idpname);
end
nexttile;
h = boxplot(intra_subjorder.SpearmanRho, intra_subjorder.method); % no default outliers
set(h,'Linew', 1.2)
hold on;
jitterAmount = 0.1;
methods      = categories(intra_subjorder.method);
for i = 1:numel(methods)
    idx = intra_subjorder.method == methods{i};
    
    x = i + (rand(sum(idx),1)-0.5)*jitterAmount; % jitter
    y = intra_subjorder.SpearmanRho(idx);
    tissues = intra_subjorder.idpname(idx);
    
    for j = 1:numel(y)
        scatter(x(j), y(j), 100, ...
            'Marker', markerMap(tissues{j}), ...
            'MarkerFaceColor', colorMap(tissues{j}), ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.8);
    end
end
hold off;
ylabel('Spearman Correlation');
ylim([-0.3,1.1]);
set(gca, 'XTickLabel', methodOrder);
title('b) Intra-scanner')
set(gca, 'FontSize', 14);
box off;
grid on;

% 2) Intra
inter_subjorder          = bigTsubjOrder(ismember(bigTsubjOrder.groupname, 'Inter-data'),:);
inter_subjorder.method   = categorical(inter_subjorder.method, methodOrder, 'Ordinal', true);
inter_subjorder.idpname  = categorical(inter_subjorder.idpname, idpOrder, 'Ordinal', true);
% Convert idpname to cell array of chars if necessary
if ~iscell(inter_subjorder.idpname)
    inter_subjorder.idpname = cellstr(inter_subjorder.idpname);
end
nexttile;
h = boxplot(inter_subjorder.SpearmanRho, inter_subjorder.method); % no default outliers
set(h,'Linew', 1.2)
hold on;
jitterAmount = 0.1;
methods      = categories(inter_subjorder.method);
for i = 1:numel(methods)
    idx = inter_subjorder.method == methods{i};
    
    x = i + (rand(sum(idx),1)-0.5)*jitterAmount; % jitter
    y = inter_subjorder.SpearmanRho(idx);
    tissues = inter_subjorder.idpname(idx);
    
    for j = 1:numel(y)
        scatter(x(j), y(j), 100, ...
            'Marker', markerMap(tissues{j}), ...
            'MarkerFaceColor', colorMap(tissues{j}), ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.8);
    end
end
hold off;
ylabel('Spearman Correlation');
ylim([-0.3,1.1]);
set(gca, 'XTickLabel', methodOrder);
title('c) Inter-scanner')
set(gca, 'FontSize', 14);
box off;
grid on;
saveas(gcf, 'withinVariability_subjorder.fig')
print('withinVariability_subjorder', '-dpng', '-r900');
close(gcf)

%% within-subject variability: RPD
rpd_struct = allResults.withinsubjvariability.rpd;
idpnames = fieldnames(rpd_struct);
allTables = cell(numel(idpnames),1);

for i = 1:numel(idpnames)
    idp = idpnames{i};
    T   = rpd_struct.(idp);
    
    % Columns to stack/melt (exclude 'groupname' and 'GroupCount')
    methodCols = setdiff(T.Properties.VariableNames, {'groupname','GroupCount'});
    
    % Stack method columns into one column 'RPDval' with method names in 'method'
    Tlong = stack(T, methodCols, ...
                  'NewDataVariableName','RPDval', ...
                  'IndexVariableName','method');
    
    % Add idpname column
    Tlong.idpname = repmat({idp}, height(Tlong), 1);
    
    allTables{i} = Tlong;
end

% Vertically concatenate all IDPs
bigRPDtable          = vertcat(allTables{:});
bigRPDtable.idpname  = cellfun(@(x) idpMap(x), bigRPDtable.idpname, 'UniformOutput', false);
bigRPDtable.method   = categorical(strrep(string(bigRPDtable.method), "_", "-"));
bigRPDtable.method   = categorical(strrep(string(bigRPDtable.method), "mean-RPD-", ""));

% 1) Repeat
repeat_rpd          = bigRPDtable(ismember(bigRPDtable.groupname, 'Repeat'),:);
repeat_rpd.method   = categorical(repeat_rpd.method, methodOrder, 'Ordinal', true);
repeat_rpd.idpname  = categorical(repeat_rpd.idpname, idpOrder, 'Ordinal', true);

% Convert idpname to cell array of chars if necessary
if ~iscell(repeat_rpd.idpname)
    repeat_rpd.idpname = cellstr(repeat_rpd.idpname);
end
f = figure;
f.Position(3:4) = [1600,1200];
tiledlayout(1,3)

nexttile;
h = boxplot(repeat_rpd.RPDval, repeat_rpd.method); % no default outliers
set(h,'Linew', 1.2)
hold on;
jitterAmount = 0.1;
methods      = categories(repeat_rpd.method);
for i = 1:numel(methods)
    idx = repeat_rpd.method == methods{i};
    
    x = i + (rand(sum(idx),1)-0.5)*jitterAmount; % jitter
    y = repeat_rpd.RPDval(idx);
    tissues = repeat_rpd.idpname(idx);
    
    for j = 1:numel(y)
        scatter(x(j), y(j), 100, ...
            'Marker', markerMap(tissues{j}), ...
            'MarkerFaceColor', colorMap(tissues{j}), ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.8);
    end
end
hold off;
ylabel('RPD (%)');
ylim([0,24]);
yticks([0:2:24])
set(gca, 'XTickLabel', methodOrder);
title('a) Repeat-scanner')
set(gca, 'FontSize', 14);
box off;
grid on;

% 2) Intra
intra_rpd          = bigRPDtable(ismember(bigRPDtable.groupname, 'Intra-scanner'),:);
intra_rpd.method   = categorical(intra_rpd.method, methodOrder, 'Ordinal', true);
intra_rpd.idpname  = categorical(intra_rpd.idpname, idpOrder, 'Ordinal', true);

% Convert idpname to cell array of chars if necessary
if ~iscell(intra_rpd.idpname)
    intra_rpd.idpname = cellstr(intra_rpd.idpname);
end
nexttile;
h = boxplot(intra_rpd.RPDval, intra_rpd.method); % no default outliers
set(h,'Linew', 1.2)
hold on;
jitterAmount = 0.1;
methods      = categories(intra_rpd.method);
for i = 1:numel(methods)
    idx = intra_rpd.method == methods{i};
    
    x = i + (rand(sum(idx),1)-0.5)*jitterAmount; % jitter
    y = intra_rpd.RPDval(idx);
    tissues = intra_rpd.idpname(idx);
    
    for j = 1:numel(y)
        scatter(x(j), y(j), 100, ...
            'Marker', markerMap(tissues{j}), ...
            'MarkerFaceColor', colorMap(tissues{j}), ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.8);
    end
end
hold off;
ylabel('RPD (%)');
ylim([0,24]);
yticks([0:2:24])
set(gca, 'XTickLabel', methodOrder);
title('b) Intra-scanner')
set(gca, 'FontSize', 14);
box off;
grid on;

% 3) Inter
inter_rpd          = bigRPDtable(ismember(bigRPDtable.groupname, 'Inter-scanner'),:);
inter_rpd.method   = categorical(inter_rpd.method, methodOrder, 'Ordinal', true);
inter_rpd.idpname  = categorical(inter_rpd.idpname, idpOrder, 'Ordinal', true);

% Convert idpname to cell array of chars if necessary
if ~iscell(inter_rpd.idpname)
    inter_rpd.idpname = cellstr(inter_rpd.idpname);
end

nexttile;
h = boxplot(inter_rpd.RPDval, inter_rpd.method); % no default outliers
set(h,'Linew', 1.2)
hold on;
jitterAmount = 0.1;
methods      = categories(inter_rpd.method);
for i = 1:numel(methods)
    idx = inter_rpd.method == methods{i};
    
    x = i + (rand(sum(idx),1)-0.5)*jitterAmount; % jitter
    y = inter_rpd.RPDval(idx);
    tissues = inter_rpd.idpname(idx);
    
    for j = 1:numel(y)
        scatter(x(j), y(j), 100, ...
            'Marker', markerMap(tissues{j}), ...
            'MarkerFaceColor', colorMap(tissues{j}), ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.8);
    end
end
hold off;
ylabel('RPD (%)');
ylim([0,24]);
yticks([0:2:24])
set(gca, 'XTickLabel', methodOrder);
title('c) Inter-scanner')
set(gca, 'FontSize', 14);
box off;
grid on;
saveas(gcf, 'withinVariability_RPD.fig')
print('withinVariability_RPD', '-dpng', '-r900');
close(gcf)