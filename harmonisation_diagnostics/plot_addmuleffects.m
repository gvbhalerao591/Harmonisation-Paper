maindir      = '/Users/psyc1586_admin/GVB_data/harmonisation_work/00_aa_skullImprovement_Oct2025/00_running_for_paper/run_eval_metrics_IQMscannerOnly/data/addmuleffects';
allItems     = dir(maindir);
dirs         = allItems([allItems.isdir]);
method_dirs  = dirs(~ismember({dirs.name}, {'.', '..'}));

%% IDP mapping
idp_list  = {'T1_SIENAX_periphGM_norm_vol',...
            'T1_SIENAX_CSF_norm_vol',...
            'T1_SIENAX_GM_norm_vol',...
            'T1_SIENAX_WM_norm_vol',...
            'T1_SIENAX_brain_norm_vol',...
            'T1_FIRST_left_hippocampus',...
            'T1_FIRST_right_hippocampus'}';
tmp_idps = {'PeriphGM', 'CSF', 'GM', 'WM', 'Brain', 'Hipp-L', 'Hipp-R'}';
mapObj   = containers.Map(idp_list, tmp_idps);

% p-value threshold
p_thr = 0.05/length(idp_list);

%% map methods

methods = {'LR_Scanner', 'LR_Site', ...
           'LME_Scanner', 'LME_ScannerFSiteR', 'LME_Site', ...
           'LME_IQM_Scanner','LME_IQM_Site', ...
           'ComBat_Scanner', 'ComBat_Site',...
           'CovBat_scanner', 'CovBat_Site', ...
           'LongComBat_Scanner', 'LongComBat_ScannerFSiteR', 'LongComBat_Site',...
           'HistMatch', ...
           'SynthSR'};

methods_new = {'LR-Scanner','LR-Site',...
                'LME-Scanner','LME-ScannerFSiteR', 'LME-Site',...
                'LME-IQM-Scanner','LME-IQM-Site', ...
                'ComBat-Scanner', 'ComBat-Site',...
                'CovBat-Scanner', 'CovBat-Site',...
                'LongComBat-Scanner', 'LongComBat-ScannerFSiteR', 'LongComBat-Site',...
                'HistMatch', ...
                'SynthSR'};              
                
% Define mapping with containers.Map
mapObj2 = containers.Map(methods, methods_new);

%%
TlongAll = table();

%%

for ii = 1:length(methods)
    method       = methods{1,ii};
    fulldir      = fullfile(maindir, method);

    add_raw      = sortrows(readtable(fullfile(fulldir,'add_raw.csv')),'Feature');
    mult_raw     = sortrows(readtable(fullfile(fulldir,'mult_raw.csv')),'Feature');
    add_harm     = sortrows(readtable(fullfile(fulldir,'add_harm.csv')),'Feature');
    mult_harm    = sortrows(readtable(fullfile(fulldir,'mult_harm.csv')),'Feature');
    
    allVals      = cell2table([add_raw.Feature, num2cell(add_raw.KRP_value), num2cell(add_harm.KRP_value),  ...
                    num2cell(mult_raw.p_value), num2cell(mult_harm.p_value)]);
    allVals.Properties.VariableNames = {'IDP', 'Add_Raw', 'Add_Harm', 'Mult_Raw', 'Mult_Harm'};
   
    % convert p-values to log scale
    tmpVals          = allVals(:,2:end);
    add_log_vals     = varfun(@(x) -log10(x), allVals, 'InputVariables', tmpVals.Properties.VariableNames);
    allVals(:,2:end) = add_log_vals;

    % Replace values in tissue column
    for i = 1:height(allVals)
        if isKey(mapObj, allVals.IDP{i})
            allVals.IDP{i} = mapObj(allVals.IDP{i});
        end
    end
    
    if sum(strcmp(methods, method)) && isKey(mapObj2, method)
        new_method = mapObj2(method);
    end

    allVals.Methodname = repmat({new_method}, height(allVals), 1);

    Tlong = stack(allVals, {'Add_Raw', 'Add_Harm', 'Mult_Raw', 'Mult_Harm'}, ...
              'NewDataVariableName','pval', ...
              'IndexVariableName','effect');
    TlongAll = [TlongAll; Tlong];
end

%% additive batch effects

add_harm = TlongAll(ismember(TlongAll.effect, 'Add_Harm'),:);
add_raw  = TlongAll(ismember(TlongAll.effect, 'Add_Raw') & ismember(TlongAll.Methodname, 'LME-Site'),:); % only raw from the "site" as batch
add_raw.Methodname(strcmp(add_raw.Methodname,'LME-Site')) = {'Raw'};
add_all  = [add_raw; add_harm];


mult_harm = TlongAll(ismember(TlongAll.effect, 'Mult_Harm'),:);
mult_raw  = TlongAll(ismember(TlongAll.effect, 'Mult_Raw') & ismember(TlongAll.Methodname, 'LME-Site'),:); % only raw from the "site" as batch
mult_raw.Methodname(strcmp(mult_raw.Methodname,'LME-Site')) = {'Raw'};
mult_all  = [mult_raw; mult_harm];


%% Define IDPs, markers, and colors
Idplist = {'Brain','CSF','GM','Hipp_L','Hipp_R','PeriphGM','WM'};
IDP_markers = {'o', '^', 's', 'v', '>', '<', 'p'};  % 7 distinct shapes
IDP_colors = [
    0.8500 0.3250 0.0980;  
    0.9290 0.6940 0.1250; 
    0.4940 0.1840 0.5560;  
    0.3010 0.7450 0.9330;  
    0.6350 0.0780 0.1840;  
    0.0000 0.4470 0.7410; 
    1.0000 0.0000 1.0000  
];

markerMap = containers.Map(Idplist, IDP_markers);
colorMap  = containers.Map(Idplist, num2cell(IDP_colors,2)); % Each row as color

%% desired order of methods and idps in plots
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
idpOrder = {'Brain', 'GM', 'PeriphGM', 'WM', 'CSF', 'Hipp_L', 'Hipp_R'};


%% Boxplot with scatter overlay

add_all.effect = [];
add_all1      = unstack(add_all, 'pval','IDP');

% Reorder rows (methods) according to methodOrder
[~, rowIdx] = ismember(methodOrder, add_all1.Methodname);
rowIdx = rowIdx(rowIdx > 0); % Keep only methods that exist
add_all1 = add_all1(rowIdx, :);

% Reorder columns (IDPs) according to idpOrder
existingIDPs = intersect(idpOrder, add_all1.Properties.VariableNames, 'stable');
add_all1 = add_all1(:, ['Methodname', existingIDPs]);

tmp_all       = table2array(add_all1(:,2:end));  % rows = methods, cols = IDPs

methodNames   = add_all1.Methodname;             % method names (first col)
idpNames      = add_all1.Properties.VariableNames(2:end); % IDPs

[numMethods, numIDPs] = size(tmp_all);

f = figure;
f.Position = [100, 100, 1200, 800];  % Large figure for A0 poster

hold on;

% Boxplot: one per method
h = boxplot(tmp_all', 'Labels', methodNames);
set(h, 'linew', 1.5);

xVals = 1:numMethods;  % x-axis positions = methods

% Store scatter handles for legend
legendHandles = gobjects(numIDPs,1);

for j = 1:numIDPs
    thisIDP = idpNames{j};
    
    % Marker & color
    if isKey(markerMap, thisIDP) && isKey(colorMap, thisIDP)
        mkr = markerMap(thisIDP);
        clr = colorMap(thisIDP);
    else
        mkr = 'o';
        clr = [0 0 0];
    end
    
    % Plot scatter for this IDP across all methods
    for i = 1:numMethods
        jitterX = (rand-0.5)*0.2; % small jitter to separate points
        h = scatter(xVals(i)+jitterX, tmp_all(i,j), 250, ...
            'Marker', mkr, ...
            'MarkerFaceColor', clr, ...
            'MarkerEdgeColor', 'k', ...
            'LineWidth', 1.5, ...
            'MarkerFaceAlpha', 0.7);
        
        % Save one handle per IDP for legend
        if i == 1
            legendHandles(j) = h;
        end
    end
end

% Add red significance threshold line
sigThresh = -log10(p_thr);
hLine = yline(sigThresh, 'r-.', 'LineWidth', 2);

hold off;
ylabel('-log10(p-value)');

% Add legend (IDPs + significance line)
legend([legendHandles; hLine], [idpNames, {sprintf('p-value=%.3f',round(p_thr,3))}], ...
       'Location', 'eastoutside','FontSize',23);

box off; grid on;
set(gca, 'GridAlpha', 0.3, 'MinorGridAlpha', 0.1);
xtickangle(45);  % Rotate x-axis labels
title('a) Additive batch effect')
set(gca, 'FontSize', 25, 'LineWidth', 1.2);
saveas(gcf, 'Figure_additive_effects.fig')
print('Figure_additive_effects', '-dpng', '-r900')
close(gcf)


%%
mult_all.effect = [];
mult_all1      = unstack(mult_all, 'pval','IDP');

% Reorder rows (methods) according to methodOrder
[~, rowIdx] = ismember(methodOrder, mult_all1.Methodname);
rowIdx = rowIdx(rowIdx > 0); % Keep only methods that exist
mult_all1 = mult_all1(rowIdx, :);

% Reorder columns (IDPs) according to idpOrder
existingIDPs = intersect(idpOrder, mult_all1.Properties.VariableNames, 'stable');
mult_all1 = mult_all1(:, ['Methodname', existingIDPs]);

tmp_all        = table2array(mult_all1(:,2:end));  % rows = methods, cols = IDPs

methodNames   = mult_all1.Methodname;             % method names (first col)
idpNames      = mult_all1.Properties.VariableNames(2:end); % IDPs

[numMethods, numIDPs] = size(tmp_all);

f = figure;
f.Position = [100, 100, 1200, 800];  % Large figure for A0 poster

hold on;

% Boxplot: one per method
h = boxplot(tmp_all', 'Labels', methodNames);
set(h, 'linew', 1.5);

xVals = 1:numMethods;  % x-axis positions = methods

% Store scatter handles for legend
legendHandles = gobjects(numIDPs,1);

for j = 1:numIDPs
    thisIDP = idpNames{j};
    
    % Marker & color
    if isKey(markerMap, thisIDP) && isKey(colorMap, thisIDP)
        mkr = markerMap(thisIDP);
        clr = colorMap(thisIDP);
    else
        mkr = 'o';
        clr = [0 0 0];
    end
    
    % Plot scatter for this IDP across all methods
    for i = 1:numMethods
        jitterX = (rand-0.5)*0.2; % small jitter to separate points
        h = scatter(xVals(i)+jitterX, tmp_all(i,j), 250, ...
            'Marker', mkr, ...
            'MarkerFaceColor', clr, ...
            'MarkerEdgeColor', 'k', ...
            'LineWidth', 1.5, ...
            'MarkerFaceAlpha', 0.7);
        
        % Save one handle per IDP for legend
        if i == 1
            legendHandles(j) = h;
        end
    end
end

% Add red significance threshold line
sigThresh = -log10(p_thr);
hLine = yline(sigThresh, 'r-.', 'LineWidth', 2);

hold off;
ylabel('-log10(p-value)');

% Add legend (IDPs + significance line)
legend([legendHandles; hLine], [idpNames, {sprintf('p-value=%.3f',round(p_thr,3))}], ...
       'Location', 'eastoutside','FontSize',23);

box off; grid on;
set(gca, 'GridAlpha', 0.3, 'MinorGridAlpha', 0.1);
xtickangle(45);  % Rotate x-axis labels
title('b) Multiplicative batch effect')
set(gca, 'FontSize', 25, 'LineWidth', 1.2);
saveas(gcf, 'Figure_multiplicative_effects.fig')
print('Figure_multiplicative_effects', '-dpng', '-r900')
close(gcf)