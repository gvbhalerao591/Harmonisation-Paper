%%

maindir = '/Users/psyc1586_admin/GVB_data/harmonisation_work/00_aa_skullImprovement_Oct2025/00_running_for_paper/run_eval_metrics_IQMscannerOnly/data';
fname   = 'data.mat';
checkid  = {"subjectID", "timepoint"}; % for sanity check

load(fullfile(maindir, fname));
allfields = fieldnames(data);

% IDP list
idp_list = {'T1_SIENAX_headsize_scaling',...
 'T1_SIENAX_periphGM_norm_vol',...
'T1_SIENAX_CSF_norm_vol',...
'T1_SIENAX_GM_norm_vol',...
'T1_SIENAX_WM_norm_vol',...
'T1_SIENAX_brain_norm_vol',...
'T1_FIRST_left_hippocampus',...
'T1_FIRST_right_hippocampus'}';

% bonferroni correction for pairwise site comparison
p_thr   = 0.05/length(idp_list);


%% Sanity check
fields   = fieldnames(data); % Get all field names

for ids = 1:length(checkid)
    sanity   = cellfun(@(f) data.(f).(checkid{1,ids}), fields, 'UniformOutput', false);% Extract all subjectID columns into a cell array
    allMatch = all(cellfun(@(x) isequal(x, sanity{1}), sanity));% Check if all are equal to the first one
    if allMatch == 1
        fprintf('\n %s match..proceeding \n', checkid{1,ids});
    else
        error('\n %s do not match..exiting \n', checkid{1,ids});
    end
end

%% Let's collate data for each IDP separately
% Configuration
refCols = {'subjectID', 'timepoint', 'Final_Age', 'Group', 'Site', 'Scanner'};
idpList = idp_list;  % List of IDPs to process

fields = fieldnames(data);
resultsStruct = struct();

% Loop through each IDP
for idpIdx = 1:length(idpList)
    idpCol = idpList{idpIdx};
    
    fprintf('\n=== Processing IDP: %s ===\n', idpCol);
    
    % Start with reference columns from Raw
    resultTable = data.Raw(:, [refCols, {idpCol}]);
    
    % Rename Raw IDP column to 'Raw'
    resultTable.Properties.VariableNames{end} = 'Raw';
    
    % Loop through other fields
    for i = 1:length(fields)
        if strcmp(fields{i}, 'Raw'), continue; end
        
        currentTable = data.(fields{i});
        colNames = currentTable.Properties.VariableNames;
        
        % Look for harmonised variants
        harmonisedPattern = ['^harmonised[A-Z]*_.*' idpCol];
        harmonisedMatch = ~cellfun(@isempty, regexpi(colNames, harmonisedPattern));
        
        if any(harmonisedMatch)
            matchedCol = colNames{find(harmonisedMatch, 1)};
            resultTable.(fields{i}) = currentTable.(matchedCol);
            fprintf('  Field "%s": Found harmonised pattern "%s"\n', fields{i}, matchedCol);
        else
            % Fallback: any column containing IDP
            idpMatch = contains(colNames, idpCol, 'IgnoreCase', true);
            if any(idpMatch)
                matchedCol = colNames{find(idpMatch, 1)};
                resultTable.(fields{i}) = currentTable.(matchedCol);
                fprintf('  Field "%s": Found IDP column "%s" (no harmonised version)\n', fields{i}, matchedCol);
            else
                warning('  Field "%s": No IDP column found!', fields{i});
            end
        end
    end
    
    % Store in results structure
    resultsStruct.(idpCol) = resultTable;
end

fprintf('\n=== Processing Complete ===\n');

%% headsize scaling - let's see what's in it

% % headsize     = resultsStruct.T1_SIENAX_headsize_scaling;
% % f = figure;
% % f.Position(3:4) = [1200, 1000];
% % tiledlayout (1,3)
% % nexttile; boxplot(headsize.Raw,      headsize.timepoint); box off; grid on; title('Raw'); ylim([0.8,1.6]), set(gca, 'FontSize', 14);
% % nexttile; boxplot(headsize.HistMatch,headsize.timepoint); box off; grid on; title('HistMatch'); ylim([0.8,1.6]), set(gca, 'FontSize', 14);
% % nexttile; boxplot(headsize.SynthSR,  headsize.timepoint); box off; grid on; title('SynthSR'); ylim([0.8,1.6]), set(gca, 'FontSize', 14);
% % saveas(gcf, 'headsize_sanity_after.fig')
% % print('headsize_sanity_after', '-dpng', 'r900')
% % close(gcf)

%% organise subjectIDs

% Sort the tables with subjectID in ascending for all fields
resultsStruct = structfun(@(t) sortrows(t, "subjectID", "ascend"), resultsStruct, 'UniformOutput', false);

% let's do sanity check again
fields   = fieldnames(resultsStruct); % Get all field names
for ids = 1:length(checkid)
    sanity   = cellfun(@(f) resultsStruct.(f).(checkid{1,ids}), fields, 'UniformOutput', false);% Extract all subjectID columns into a cell array
    allMatch = all(cellfun(@(x) isequal(x, sanity{1}), sanity));% Check if all are equal to the first one
    if allMatch == 1
        fprintf('\n %s match..proceeding \n', checkid{1,ids});
    else
        error('\n %s do not match..exiting \n', checkid{1,ids});
    end
end

%%

finalidplist = {...
 'T1_SIENAX_periphGM_norm_vol',...
'T1_SIENAX_CSF_norm_vol',...
'T1_SIENAX_GM_norm_vol',...
'T1_SIENAX_WM_norm_vol',...
'T1_SIENAX_brain_norm_vol',...
'T1_FIRST_left_hippocampus',...
'T1_FIRST_right_hippocampus'}';

methods = fieldnames(data);

subjectID = "subjectID";
Group     = "Group";
timepoint = "timepoint";
age       = "Final_Age";
site      = "Site";
scanner   = "Scanner";

%% 1) RPD - within-subject variability

for idp = 1:length(finalidplist)
    
    tmpidp   = finalidplist{idp,1};
    all_data = resultsStruct.(tmpidp);
    
    rpdFun   = @(x) abs(x(1) - x(2)) / mean(x) * 100; % ratio of absolute difference by mean 
    
    % Apply it using varfun
    RPDTable = varfun(rpdFun, all_data, ...
        'InputVariables', methods, ...
        'GroupingVariables', 'subjectID');

    % Rename the result column for clarity
    RPDTable.Properties.VariableNames = strrep(RPDTable.Properties.VariableNames,'Fun_','RPD_');
    
    uniqueRows         = unique(all_data(:, [subjectID, Group]), 'rows', 'stable');
    RPDTable.groupname = uniqueRows.(Group);
    
    % sort table
    RPDTable.GroupType = categorical(RPDTable.groupname, ...
        {'Repeat', 'Intra-scanner', 'Inter-scanner'}, ...
        'Ordinal', true);
    
    RPDTable_sorted    = sortrows(RPDTable, 'GroupType');
    allRPD.(tmpidp)    = RPDTable_sorted;
    meanRPD            = varfun(@mean, RPDTable_sorted, ...
                                'InputVariables', RPDTable_sorted.Properties.VariableNames(3:end-2), ...
                                'GroupingVariables', 'groupname');
    % meanRPD across subject - all
    tmp_all      = table2array(varfun(@mean, RPDTable_sorted,'InputVariables',RPDTable_sorted.Properties.VariableNames(3:end-2)));
    all_data_tmp = cell2table(['All-data', num2cell(sum(meanRPD.GroupCount)), num2cell(tmp_all)],'VariableNames',meanRPD.Properties.VariableNames);
  
    % Append
    avgRPD.(tmpidp) = [meanRPD; all_data_tmp];
end

%% subject order consistency - within-subject variability
nPerm = 10000;

for idp = 1:length(finalidplist)
    
    tmpidp   = finalidplist{idp,1};
    all_data = resultsStruct.(tmpidp);
    
    tp0_all_data    = all_data(ismember(all_data.(timepoint),'TP0'),:);
    tp1_all_data    = all_data(ismember(all_data.(timepoint),'TP1'),:);
    subjOrderTable0 = evaluate_subjOrder(tp0_all_data, tp1_all_data, methods, nPerm);
    subjOrderTable0 = [subjOrderTable0(:,1), table(repmat("All-data", height(subjOrderTable0), 1), 'VariableNames', {'groupname'}), subjOrderTable0(:,2:end)];


    tp0_all_data    = all_data(ismember(all_data.(timepoint),'TP0') & ismember(all_data.(Group),'Inter-scanner'),:);
    tp1_all_data    = all_data(ismember(all_data.(timepoint),'TP1') & ismember(all_data.(Group),'Inter-scanner'),:);
    subjOrderTable1 = evaluate_subjOrder(tp0_all_data, tp1_all_data, methods, nPerm);
    subjOrderTable1 = [subjOrderTable1(:,1), table(repmat("Inter-data", height(subjOrderTable1), 1), 'VariableNames', {'groupname'}), subjOrderTable1(:,2:end)];

    tp0_all_data    = all_data(ismember(all_data.(timepoint),'TP0') & ismember(all_data.(Group),'Intra-scanner'),:);
    tp1_all_data    = all_data(ismember(all_data.(timepoint),'TP1') & ismember(all_data.(Group),'Intra-scanner'),:);
    subjOrderTable2 = evaluate_subjOrder(tp0_all_data, tp1_all_data, methods, nPerm);
    subjOrderTable2 = [subjOrderTable2(:,1), table(repmat("Intra-data", height(subjOrderTable2), 1), 'VariableNames', {'groupname'}), subjOrderTable2(:,2:end)];

    tp0_all_data    = all_data(ismember(all_data.(timepoint),'TP0') & ismember(all_data.(Group),'Repeat'),:);
    tp1_all_data    = all_data(ismember(all_data.(timepoint),'TP1') & ismember(all_data.(Group),'Repeat'),:);
    subjOrderTable3 = evaluate_subjOrder(tp0_all_data, tp1_all_data, methods, nPerm);
    subjOrderTable3 = [subjOrderTable3(:,1), table(repmat("Repeat", height(subjOrderTable3), 1), 'VariableNames', {'groupname'}), subjOrderTable3(:,2:end)];

    subjOrder_tables.(tmpidp)  = [subjOrderTable0; subjOrderTable1; subjOrderTable2; subjOrderTable3];
end

%% models on raw and harmonised data
% 1) Get batch effect and pairwise significance
% 2) Get ICC - between subject variability
% 3) Get association with age and timepoint

for idp = 1:length(finalidplist)
    
    tmpidp   = finalidplist{idp,1};
    all_data = resultsStruct.(tmpidp);

    all_data.age = all_data.(age);
    
    % make timepoint also categorical
    all_data.timepoint = categorical(all_data.(timepoint));
    all_data.scanner   = categorical(all_data.(scanner));
    all_data.site      = categorical(all_data.(site));

    % zscore volumes and age 
    tmp_methods = ['age',methods'];
    zscore_vols = varfun(@zscore,all_data,"InputVariables",tmp_methods);
    all_data    = [all_data, zscore_vols];

    % Get the count of samples per site
    site_counts = countcats(all_data.site);
    site_names  = categories(all_data.site);
    
    % Find the site with the most samples
    [~, idx_max] = max(site_counts);
    ref_site     = site_names{idx_max};
    
    % Reorder categories to make the most frequent site the reference
    new_order     = [ref_site; setdiff(site_names, ref_site, 'stable')];
    all_data.site = reordercats(all_data.site, new_order); 

    %% Fit diagnostic models - site
    
    outs_mdl_site = cell(length(methods), 18);

    for ii = 1:length(methods)
        
        responseVar = ['zscore_' methods{ii,1}];

        ageVar      = 'zscore_age';
        timeVar     = 'timepoint';
        batchVar    = 'site';
        ranEff      = '(1|subjectID)';
        
        outs_mdl_site{ii,1} = strrep(methods{ii,1},'_','-');
        outs_mdl_site{ii,2} = 'site';
        outs_mdl_site{ii,3} = tmpidp;

        % To see residual batch effects - site
        formula_str          = sprintf('%s ~ %s + %s + %s + %s',responseVar, ageVar, timeVar, batchVar, ranEff);
        mdl1                 = fitlme(all_data, formula_str);
        outs_mdl_site{ii,4}  = pairwise_siteTests(mdl1, 'site', all_data,p_thr);
        tmp                  = anova(mdl1);
        outs_mdl_site{ii,5}  = sum(double(tmp(3, 5))<0.05);

        % To get ICC
        formula_str           = sprintf('%s ~ 1 + %s', responseVar, ranEff);
        mdl2                  = fitlme(all_data, formula_str);
        tmp_subj_variance     =  cell2mat(mdl2.covarianceParameters);   
        outs_mdl_site{ii,6}   = 0;
        outs_mdl_site{ii,7}   = tmp_subj_variance;
        residual_variance     = mdl2.MSE;
        outs_mdl_site{ii,8}   = residual_variance; 
        outs_mdl_site{ii,9}   = tmp_subj_variance / (tmp_subj_variance + residual_variance);
        outs_mdl_site{ii,10}  = residual_variance/tmp_subj_variance;

        % To get age and time effects
        formula_str            = sprintf('%s ~ %s + %s + %s', responseVar, ageVar, timeVar, ranEff);
        mdl3                   = fitlme(all_data, formula_str);
        coeffs                 = mdl3.Coefficients;
        outs_mdl_site{ii,11}   = coeffs.Estimate(strcmp(coeffs.Name, 'zscore_age'));
        outs_mdl_site{ii,12}   = coeffs.pValue(strcmp(coeffs.Name,   'zscore_age'));
        outs_mdl_site{ii,13}   = coeffs.Lower(strcmp(coeffs.Name,    'zscore_age'));
        outs_mdl_site{ii,14}   = coeffs.Upper(strcmp(coeffs.Name,    'zscore_age'));
        outs_mdl_site{ii,15}   = coeffs.Estimate(strcmp(coeffs.Name, 'timepoint_TP1'));
        outs_mdl_site{ii,16}   = coeffs.pValue(strcmp(coeffs.Name,   'timepoint_TP1'));
        outs_mdl_site{ii,17}   = coeffs.Lower(strcmp(coeffs.Name,    'timepoint_TP1'));
        outs_mdl_site{ii,18}   = coeffs.Upper(strcmp(coeffs.Name,    'timepoint_TP1'));
    end

    outs_mdl_site_tab = cell2table(outs_mdl_site, 'VariableNames',{'method','Effects','idpname','n_is_batchSig','anova_batches','Site_Var','Subj_Var','Resid_Var','ICC','WCV','est_Age','pval_Age','ciL_Age','ciU_Age','est_time','pval_time','ciL_time','ciU_time'});
    
    %% Fit diagnostic models - scanner

    outs_mdl_scanner = cell(length(methods), 18);

    for ii = 1:length(methods)

        responseVar = ['zscore_' methods{ii,1}];
         
        ageVar      = 'zscore_age';
        timeVar     = 'timepoint';
        batchVar    = 'scanner';
        ranEff      = '(1|subjectID)';
        
        outs_mdl_scanner{ii,1} = strrep(methods{ii,1},'_','-');
        outs_mdl_scanner{ii,2} = 'scanner';
        outs_mdl_scanner{ii,3} = tmpidp;

        % To see residual batch effects - scanner
        formula_str             = sprintf('%s ~ %s + %s + %s + %s',responseVar, ageVar, timeVar, batchVar, ranEff);
        mdl1                    = fitlme(all_data, formula_str);
        outs_mdl_scanner{ii,4}  = pairwise_siteTests(mdl1, 'scanner', all_data, p_thr);
        tmp                     = anova(mdl1);
        outs_mdl_scanner{ii,5}  = sum(double(tmp(3, 5))<0.05);

        % To get ICC
        formula_str                   = sprintf('%s ~ 1 + %s', responseVar, ranEff);
        mdl2                          = fitlme(all_data, formula_str);
        tmp_subj_variance             =  cell2mat(mdl2.covarianceParameters);   
        outs_mdl_scanner{ii,6}        = 0;
        outs_mdl_scanner{ii,7}        = tmp_subj_variance;
        residual_variance             = mdl2.MSE;
        outs_mdl_scanner{ii,8}        = residual_variance; 
        outs_mdl_scanner{ii,9}        = tmp_subj_variance / (tmp_subj_variance + residual_variance);
        outs_mdl_scanner{ii,10}       = residual_variance/tmp_subj_variance;

        % To get age and time effects
        formula_str                = sprintf('%s ~ %s + %s + %s', responseVar, ageVar, timeVar, ranEff);
        mdl3                       = fitlme(all_data, formula_str);
        coeffs                     = mdl3.Coefficients;
        outs_mdl_scanner{ii,11}    = coeffs.Estimate(strcmp(coeffs.Name, 'zscore_age'));
        outs_mdl_scanner{ii,12}    = coeffs.pValue(strcmp(coeffs.Name,   'zscore_age'));
        outs_mdl_scanner{ii,13}    = coeffs.Lower(strcmp(coeffs.Name,    'zscore_age'));
        outs_mdl_scanner{ii,14}    = coeffs.Upper(strcmp(coeffs.Name,    'zscore_age'));
        outs_mdl_scanner{ii,15}    = coeffs.Estimate(strcmp(coeffs.Name, 'timepoint_TP1'));
        outs_mdl_scanner{ii,16}    = coeffs.pValue(strcmp(coeffs.Name,   'timepoint_TP1'));
        outs_mdl_scanner{ii,17}    = coeffs.Lower(strcmp(coeffs.Name,    'timepoint_TP1'));
        outs_mdl_scanner{ii,18}    = coeffs.Upper(strcmp(coeffs.Name,    'timepoint_TP1'));
    end

    outs_mdl_scanner_tab = cell2table(outs_mdl_scanner, 'VariableNames',{'method','Effects','idpname','n_is_batchSig','anova_batches','Site_Var','Subj_Var','Resid_Var','ICC','WCV','est_Age','pval_Age','ciL_Age','ciU_Age','est_time','pval_time','ciL_time','ciU_time'});
    
    outs_mdl.(tmpidp) = [outs_mdl_site_tab; outs_mdl_scanner_tab];
end

bigTable = [];
for i = 1:length(finalidplist)
    bigTable = [bigTable; outs_mdl.(finalidplist{i})];
end

common_cols              = {'method','Effects','idpname'};
residualVariability_cols = bigTable(:,[common_cols, {'n_is_batchSig','anova_batches'}]);
betSubjVariability_cols  = bigTable(:,[common_cols, {'Subj_Var','Resid_Var','ICC','WCV'}]);
bioAgeVariability_cols   = bigTable(:,[common_cols, {'est_Age','pval_Age','ciL_Age','ciU_Age'}]);
bioTimeVariability_cols  = bigTable(:,[common_cols, {'est_time','pval_time','ciL_time','ciU_time'}]);

%% Multivariate - mahalanobis distance
mdidplist = {...
 'T1_SIENAX_periphGM_norm_vol',...
'T1_SIENAX_CSF_norm_vol',...
'T1_SIENAX_GM_norm_vol',...
'T1_SIENAX_WM_norm_vol',...
'T1_SIENAX_brain_norm_vol',...
'T1_FIRST_left_hippocampus',...
'T1_FIRST_right_hippocampus'}';

methodfields = fieldnames(data);

% We need reference for site/scanner info so first field is reference which is Raw here
tmpData           = data.(methodfields{1,1});
siteVar_double    = double(categorical(tmpData.(site)));
ourScanners       = double(categorical(tmpData.(scanner)));
MD_outs           = zeros(length(unique(siteVar_double)), length(methodfields));


for mm = 1:length(methodfields)
    
    tmpfield = methodfields{mm,1};
    dataTab  = data.(tmpfield);

    dataTab_forMD = dataTab(:, contains(dataTab.Properties.VariableNames, 'harmonised'));

    if isempty(dataTab_forMD)
       dataTab_forMD = dataTab(:, contains(dataTab.Properties.VariableNames, mdidplist));
    end

    MD_outs(:,mm)  = get_MD_numerically_stable(dataTab_forMD,  siteVar_double);
end
MD_outs(length(unique(siteVar_double))+1,:) = mean(MD_outs);
siteVar_cat        = cellstr(unique(categorical(tmpData.(site))));
siteVar_cat(end+1) = {'average_batch'};
all_MDs            = cell2table([siteVar_cat, num2cell(MD_outs)],'VariableNames',['batch',methodfields']);

%% Save results

allResults.input.idpwise                   = resultsStruct;
allResults.input.data_methodwise           = data;
allResults.withinsubjvariability.rpd       = avgRPD;
allResults.withinsubjvariability.subjorder = subjOrder_tables;
allResults.residualvariability             = residualVariability_cols;
allResults.betSubjVariability              = betSubjVariability_cols;
allResults.bioAgeVariability               = bioAgeVariability_cols;
allResults.bioTimeVariability              = bioTimeVariability_cols;
allResults.multivariateVariability         = all_MDs;

save(fullfile(maindir, 'results.mat'), 'allResults')





