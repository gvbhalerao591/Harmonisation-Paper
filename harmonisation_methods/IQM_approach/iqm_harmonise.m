function data = iqm_harmonise_working(data, idp_list, qc_list, iqm_variance, p_thr, max_qcs,outfilename)
    % iqm_harmonise_working: Harmonization matching the proven working code
    %
    % KEY DIFFERENCES from "rigorous" versions:
    % 1. Uses ORIGINAL volume (not additive-corrected) for multiplicative QC selection
    % 2. Multiplicative criteria: QC should NOT be strongly variance-driven (opposite!)
    % 3. Always attempts additive correction (no pre-gating)
    % 4. Model comparison ONLY for multiplicative
    %
    % Inputs:
    %   data         - Table containing volumes, covariates, and QC metrics
    %   idp_list     - Cell array of volume names to harmonize
    %   qc_list      - Cell array of QC metric names to use
    %   iqm_variance - Variance threshold for PCA (e.g., 95)
    %   p_thr        - P-value threshold (e.g., 0.05)
    %   max_qcs      - Maximum QCs per volume (default: Inf, no limit)
    %
    % Outputs:
    %   data         - Table with harmonized volumes added
    
    if nargin < 6
        max_qcs = Inf;
    end
    
    fprintf('=== IQM Harmonization (Working Version) ===\n');
    fprintf('P-value threshold: %.4f\n', p_thr);
    fprintf('Max QCs per volume: %s\n', num2str(max_qcs));
    fprintf('PCA variance threshold: %.1f%%\n\n', iqm_variance);
    
    %% ========================================================================
    %% STEP 1: PCA ON QC METRICS
    %% ========================================================================
    
    fprintf('Step 1: Preparing QC metrics and PCA...\n');
    
    if ~iscell(qc_list)
        error('qc_list must be a cell array');
    end
    
    % Validate and extract QCs
    missing_qcs = ~ismember(qc_list, data.Properties.VariableNames);
    if any(missing_qcs)
        error('Missing QC metrics: %s', strjoin(qc_list(missing_qcs), ', '));
    end
    
    data_iqms = data(:, qc_list);
    tmp_qc = table2array(data_iqms);
    std_tmp_qc = std(tmp_qc)';
    remove_qc = find(std_tmp_qc == 0);
    
    if ~isempty(remove_qc)
        fprintf('  Removing %d zero-variance QCs\n', length(remove_qc));
        data_iqms(:, remove_qc) = [];
        qc_list(remove_qc) = [];
    end
    
    % PCA
    qc_matrix = table2array(data_iqms);
    qc_z = zscore(qc_matrix);
    [coeff_iqm, scores_iqm, latent_iqm, ~, explained_iqm] = pca(qc_z);
    cumulative_explained_iqm = cumsum(explained_iqm);
    num_components_iqm = find(cumulative_explained_iqm >= iqm_variance, 1);
    
    fprintf('  Retained %d PCs (%.1f%% variance)\n', ...
            num_components_iqm, cumulative_explained_iqm(num_components_iqm));
    
    pca_names = arrayfun(@(i) sprintf('QC%d', i), 1:num_components_iqm, 'UniformOutput', false);
    qc_scores = array2table(scores_iqm(:, 1:num_components_iqm), 'VariableNames', pca_names);
    
    % Merge data
    data = [data, qc_scores];
    data.batch = categorical(data.batch);
    data.timepoint = categorical(data.timepoint);
    data.subjectID = categorical(data.subjectID);
    
    total_pcs = size(qc_scores, 2);
    
    %% ========================================================================
    %% STEP 2: ADDITIVE QC SELECTION
    %% ========================================================================
    
    fprintf('\nStep 2: Selecting QCs for ADDITIVE correction...\n');
    fprintf('  Criteria: batch-driven, not age/volume-driven\n');
    
    all_out_add = cell(length(idp_list), 2);
    
    for vols = 1:length(idp_list)
        volname = idp_list{vols};
        outs = cell(total_pcs, 6);
        
        for pcs = 1:total_pcs
            qcname = sprintf('QC%d', pcs);
            formula1 = sprintf('%s ~ age + timepoint + batch + %s + (1|subjectID)', qcname, volname);
            formula2 = sprintf('age ~ %s + timepoint + batch + %s + (1|subjectID)', qcname, volname);
            
            try
                mdl_qc = fitlme(data, formula1);
                mdl_age = fitlme(data, formula2);
                tmp_mdl_qc = anova(mdl_qc);
                tmp_mdl_age = anova(mdl_age);
                
                age_idx = find(strcmp(tmp_mdl_qc.Term, 'age'));
                vol_idx = find(strcmp(tmp_mdl_qc.Term, volname));
                batch_idx = find(strcmp(tmp_mdl_qc.Term, 'batch'));
                qc_idx = find(strcmp(tmp_mdl_age.Term, qcname));
                
                pval_age = tmp_mdl_qc.pValue(age_idx);
                pval_vol = tmp_mdl_qc.pValue(vol_idx);
                pval_batch = tmp_mdl_qc.pValue(batch_idx);
                pval_qc = tmp_mdl_age.pValue(qc_idx);
                
                pass_filter = (pval_age > p_thr) && (pval_vol > p_thr) && ...
                             (pval_batch < p_thr) && (pval_qc > p_thr);
                
                outs{pcs, 1} = qcname;
                outs{pcs, 2} = pval_age;
                outs{pcs, 3} = pval_vol;
                outs{pcs, 4} = pval_batch;
                outs{pcs, 5} = pval_qc;
                outs{pcs, 6} = pass_filter;
            catch ME
                warning('Error in additive selection for %s, %s: %s', volname, qcname, ME.message);
                outs{pcs, 1} = qcname;
                [outs{pcs, 2:6}] = deal(NaN, NaN, NaN, NaN, 0);
            end
        end
        
        all_out_add{vols, 1} = volname;
        all_out_add{vols, 2} = outs;
    end
    
    % Build additive mask
    good_pcs_add = zeros(total_pcs, length(idp_list));
    for idps = 1:length(idp_list)
        tmp_good = cell2table(all_out_add{idps, 2});
        good_pcs_add(:, idps) = table2array(tmp_good(:, end));
    end
    
    n_qcs_add = sum(good_pcs_add, 1);
    fprintf('  QCs selected: min=%d, max=%d, mean=%.1f\n', ...
            min(n_qcs_add), max(n_qcs_add), mean(n_qcs_add));
    
    %% ========================================================================
    %% STEP 3: MULTIPLICATIVE QC SELECTION (FROM ORIGINAL VOLUMES!)
    %% ========================================================================
    
    fprintf('\nStep 3: Selecting QCs for MULTIPLICATIVE correction...\n');
    fprintf('  KEY: Using ORIGINAL volumes (not additive-corrected) for variance proxy\n');
    fprintf('  Criteria: batch-driven, NOT strongly variance-driven\n');
    
    all_out_mult = cell(length(idp_list), 2);
    
    for vols = 1:length(idp_list)
        volname = idp_list{vols};
        outs = cell(total_pcs, 6);
        
        % CRITICAL: Use ORIGINAL volume for biology-only residuals
        lme_bio = fitlme(data, sprintf('%s ~ age + timepoint + (1|subjectID)', volname));
        res = residuals(lme_bio);
        data.log_r2_tmp = log(max(res.^2, eps));
        % data.log_r2_tmp = log(data.(volname));
        
        for pcs = 1:total_pcs
            qcname = sprintf('QC%d', pcs);
            
            try
                % Variance model: QC ~ age + timepoint + batch + log_r2
                Tvar = data(:, {'age', 'timepoint', 'batch', 'subjectID', qcname, 'log_r2_tmp'});
                mdl_var = fitlme(Tvar, sprintf('%s ~ age + timepoint + batch + log_r2_tmp + (1|subjectID)', qcname));
                tmp_mdl_var = anova(mdl_var);
                
                age_idx = find(strcmp(tmp_mdl_var.Term, 'age'));
                vol_idx = find(strcmp(tmp_mdl_var.Term, 'log_r2_tmp'));
                batch_idx = find(strcmp(tmp_mdl_var.Term, 'batch'));
                
                pval_age = tmp_mdl_var.pValue(age_idx);
                pval_vol_var = tmp_mdl_var.pValue(vol_idx);
                pval_batch = tmp_mdl_var.pValue(batch_idx);
                
                % Age guard
                mdl_age = fitlme(data, sprintf('age ~ %s + timepoint + batch + log_r2_tmp + (1|subjectID)', qcname));
                tmp_mdl_age = anova(mdl_age);
                qc_idx = find(strcmp(tmp_mdl_age.Term, qcname));
                pval_qc = tmp_mdl_age.pValue(qc_idx);
                
                % CRITICAL DIFFERENCE: pval_vol_var > p_thr (NOT <)
                % QC should capture batch but NOT be variance-driven
                pass_filter = (pval_age > p_thr) && (pval_vol_var > p_thr) && ...
                             (pval_batch < p_thr) && (pval_qc > p_thr);
                
                outs{pcs, 1} = qcname;
                outs{pcs, 2} = pval_age;
                outs{pcs, 3} = pval_vol_var;
                outs{pcs, 4} = pval_batch;
                outs{pcs, 5} = pval_qc;
                outs{pcs, 6} = pass_filter;
            catch ME
                warning('Error in multiplicative selection for %s, %s: %s', volname, qcname, ME.message);
                outs{pcs, 1} = qcname;
                [outs{pcs, 2:6}] = deal(NaN, NaN, NaN, NaN, 0);
            end
        end
        
        all_out_mult{vols, 1} = volname;
        all_out_mult{vols, 2} = outs;
        data.log_r2_tmp = [];
    end
    
    % Build multiplicative mask
    good_pcs_mult = zeros(total_pcs, length(idp_list));
    for idps = 1:length(idp_list)
        tmp_good = cell2table(all_out_mult{idps, 2});
        good_pcs_mult(:, idps) = table2array(tmp_good(:, end));
    end
    
    n_qcs_mult = sum(good_pcs_mult, 1);
    fprintf('  QCs selected: min=%d, max=%d, mean=%.1f\n', ...
            min(n_qcs_mult), max(n_qcs_mult), mean(n_qcs_mult));
    
    %% ========================================================================
    %% STEP 4: APPLY CORRECTIONS (SEQUENTIAL: ADDITIVE → MULTIPLICATIVE)
    %% ========================================================================
    
    fprintf('\nStep 4: Applying corrections...\n');
    
    qc_selection = struct();
    qc_selection.volumes = idp_list;
    qc_selection.additive_qcs = cell(length(idp_list), 1);
    qc_selection.multiplicative_qcs = cell(length(idp_list), 1);
    qc_selection.additive_count = zeros(length(idp_list), 1);
    qc_selection.multiplicative_count = zeros(length(idp_list), 1);
    qc_selection.multiplicative_model_pval = nan(length(idp_list), 1);
    qc_selection.multiplicative_applied = false(length(idp_list), 1);
    
    y_harmonised = zeros(size(data, 1), length(idp_list));
    
    for v = 1:length(idp_list)
        volname = idp_list{v};
        
        %% ADDITIVE CORRECTION (always attempted)
        qc_mask_add = find(good_pcs_add(:, v) == 1);
        
        % Apply max_qcs limit
        if ~isinf(max_qcs) && length(qc_mask_add) > max_qcs
            qc_mask_add = qc_mask_add(1:max_qcs);
        end
        
        if isempty(qc_mask_add)
            % No additive correction
            data.tmp_noAdditive = data.(volname);
            qc_selection.additive_qcs{v} = {};
            qc_selection.additive_count(v) = 0;
        else
            % Apply additive correction
            qc_vars_add = arrayfun(@(i) sprintf('QC%d', i), qc_mask_add, 'UniformOutput', false);
            qc_terms_add = strjoin(qc_vars_add, ' + ');
            
            formula_additive = sprintf('%s ~ age + timepoint + %s + (1|subjectID)', volname, qc_terms_add);
            lme_additive = fitlme(data, formula_additive);
            
            beta_add = fixedEffects(lme_additive);
            beta_names_add = lme_additive.CoefficientNames;
            idx_qc_add = ismember(beta_names_add, qc_vars_add);
            X_qc_add = data{:, beta_names_add(idx_qc_add)};
            additive_effect = X_qc_add * beta_add(idx_qc_add);
            
            data.tmp_noAdditive = data.(volname) - additive_effect;
            qc_selection.additive_qcs{v} = qc_vars_add;
            qc_selection.additive_count(v) = length(qc_mask_add);
        end
        
        % Ensure positive values
        epsilon = 1e-6;
        data.tmp_noAdditive = max(data.tmp_noAdditive, epsilon);
        
        %% MULTIPLICATIVE CORRECTION (with model comparison)
        qc_mask_mult = find(good_pcs_mult(:, v) == 1);
        
        % Apply max_qcs limit
        if ~isinf(max_qcs) && length(qc_mask_mult) > max_qcs
            qc_mask_mult = qc_mask_mult(1:max_qcs);
        end
        
        if isempty(qc_mask_mult)
            % No multiplicative correction
            y_harmonised(:, v) = data.tmp_noAdditive;
            qc_selection.multiplicative_qcs{v} = {};
            qc_selection.multiplicative_count(v) = 0;
            qc_selection.multiplicative_applied(v) = false;
        else
            % Test if QCs improve model
            qc_vars_mult = arrayfun(@(i) sprintf('QC%d', i), qc_mask_mult, 'UniformOutput', false);
            qc_terms_mult = strjoin(qc_vars_mult, ' + ');
            
            data.log_tmp_noAdditive = log(data.tmp_noAdditive);
            
            formula_reduced = sprintf('log_tmp_noAdditive ~ age + timepoint + (1|subjectID)');
            formula_full = sprintf('log_tmp_noAdditive ~ age + timepoint + %s + (1|subjectID)', qc_terms_mult);
            
            try
                lme_reduced = fitlme(data, formula_reduced);
                lme_full = fitlme(data, formula_full);
                cmp = compare(lme_reduced, lme_full);
                model_pval = cmp.pValue(2);
                
                qc_selection.multiplicative_model_pval(v) = model_pval;
                
                if model_pval < p_thr
                    % Apply multiplicative correction
                    beta_mult = fixedEffects(lme_full);
                    beta_names_mult = lme_full.CoefficientNames;
                    idx_qc_mult = ismember(beta_names_mult, qc_vars_mult);
                    X_qc_mult = data{:, beta_names_mult(idx_qc_mult)};
                    
                    lp = X_qc_mult * beta_mult(idx_qc_mult);
                    lp_centered = lp - mean(lp, 'omitnan');
                    multiplicative_effect = exp(lp_centered);
                    
                    y_harmonised(:, v) = data.tmp_noAdditive ./ multiplicative_effect;
                    qc_selection.multiplicative_qcs{v} = qc_vars_mult;
                    qc_selection.multiplicative_count(v) = length(qc_mask_mult);
                    qc_selection.multiplicative_applied(v) = true;
                    
                    fprintf('  %s: Add=%d QCs, Mult=%d QCs (p=%.4f)\n', ...
                            volname, length(qc_mask_add), length(qc_mask_mult), model_pval);
                else
                    % QCs don't improve model - skip multiplicative
                    y_harmonised(:, v) = data.tmp_noAdditive;
                    qc_selection.multiplicative_qcs{v} = {};
                    qc_selection.multiplicative_count(v) = 0;
                    qc_selection.multiplicative_applied(v) = false;
                end
            catch ME
                warning('Model comparison failed for %s: %s', volname, ME.message);
                y_harmonised(:, v) = data.tmp_noAdditive;
                qc_selection.multiplicative_qcs{v} = {};
                qc_selection.multiplicative_count(v) = 0;
                qc_selection.multiplicative_applied(v) = false;
            end
            
            data.log_tmp_noAdditive = [];
        end
        
        % Store harmonised result
        harmonised_colname = ['harmonised_', volname];
        data.(harmonised_colname) = y_harmonised(:, v);
    end
    
    % Cleanup
    if ismember('tmp_noAdditive', data.Properties.VariableNames)
        data.tmp_noAdditive = [];
    end
    
    %% ========================================================================
    %% STEP 5: SAVE AND SUMMARIZE
    %% ========================================================================
    
    fprintf('\nStep 5: Saving results...\n');
    
    writetable(data, outfilename);
    fprintf('  Data saved: iqm_harmonized_working.csv\n');
    
    qc_selection.p_threshold = p_thr;
    qc_selection.max_qcs = max_qcs;
    qc_selection.pca_info = struct();
    qc_selection.pca_info.num_components = num_components_iqm;
    qc_selection.pca_info.variance_explained = cumulative_explained_iqm(num_components_iqm);
    
    save('qc_selection_working.mat', 'qc_selection', 'all_out_add', 'all_out_mult', ...
         'good_pcs_add', 'good_pcs_mult');
    fprintf('  QC selection saved: qc_selection_working.mat\n');
    
    %% Summary
    fprintf('\n=== SUMMARY ===\n');
    fprintf('Additive correction:\n');
    fprintf('  Applied to: %d/%d volumes\n', sum(qc_selection.additive_count > 0), length(idp_list));
    if sum(qc_selection.additive_count > 0) > 0
        fprintf('  QCs used: min=%d, max=%d, mean=%.1f\n', ...
                min(qc_selection.additive_count(qc_selection.additive_count > 0)), ...
                max(qc_selection.additive_count), ...
                mean(qc_selection.additive_count(qc_selection.additive_count > 0)));
    end
    
    fprintf('\nMultiplicative correction:\n');
    fprintf('  QCs selected for: %d/%d volumes\n', sum(qc_selection.multiplicative_count > 0), length(idp_list));
    fprintf('  Actually applied: %d/%d volumes\n', sum(qc_selection.multiplicative_applied), length(idp_list));
    if sum(qc_selection.multiplicative_applied) > 0
        applied_idx = qc_selection.multiplicative_applied;
        fprintf('  QCs used: min=%d, max=%d, mean=%.1f\n', ...
                min(qc_selection.multiplicative_count(applied_idx)), ...
                max(qc_selection.multiplicative_count(applied_idx)), ...
                mean(qc_selection.multiplicative_count(applied_idx)));
        fprintf('  Model p-values: min=%.4f, max=%.4f, median=%.4f\n', ...
                min(qc_selection.multiplicative_model_pval(applied_idx)), ...
                max(qc_selection.multiplicative_model_pval(applied_idx)), ...
                median(qc_selection.multiplicative_model_pval(applied_idx)));
    end
    
    fprintf('\n=== KEY FEATURES (MATCHING WORKING CODE) ===\n');
    fprintf('✓ Multiplicative QCs selected from ORIGINAL volumes\n');
    fprintf('✓ Multiplicative criteria: batch-driven, NOT variance-driven\n');
    fprintf('✓ Always attempts additive (no pre-gating)\n');
    fprintf('✓ Model comparison only for multiplicative\n');
    fprintf('\nDone!\n');
    
end