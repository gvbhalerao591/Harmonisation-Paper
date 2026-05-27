function [data, qc_selection, add_detail_tbl, mult_detail_tbl, summary_tbl] = iqm_harmonise_flexible( ...
    data, idp_list, qc_list, preserve_covars, adjust_covars, reverse_guard_covars, categorical_covars, ...
    batch_col, subject_col, age_source_col, batch_source_col, age_col, iqm_variance, p_thr, max_qcs, apply_pca, outfilename)
% IQM harmonisation with flexible covariates and optional PCA/raw QC basis
%
% This version is intentionally flat:
%   - no local helper functions
%   - additive and multiplicative screening use fitlme + anova
%   - compare(...) is used only for the final multiplicative apply/skip gate
%
% Recommended use:
%   - pass empty lists as {}
%   - use {} for adjust_covars when you want none
%
% Example:
%   [data_out, qc_selection, add_detail_tbl, mult_detail_tbl, summary_tbl] = iqm_harmonise_flexible( ...
%       data, idp_list, iqm_list, {'age','sex'}, {}, {'age'}, {'scan_session','sex'}, ...
%       'scan_session', 'subject', 'age', 'scan_session', 'age', 95, 0.05, Inf, true, 'iqm_harmonised.csv');

    %% ---------------------------------------------------------------------
    %% Defaults
    %% ---------------------------------------------------------------------
    if nargin < 4 || isempty(preserve_covars), preserve_covars = {'age'}; end
    if nargin < 5
        adjust_covars = {'timepoint'};
    end
    if nargin < 6, reverse_guard_covars = {}; end
    if nargin < 7 || isempty(categorical_covars), categorical_covars = {'timepoint', 'batch'}; end
    if nargin < 8 || isempty(batch_col), batch_col = 'batch'; end
    if nargin < 9 || isempty(subject_col), subject_col = 'subjectID'; end
    if nargin < 10 || isempty(age_source_col), age_source_col = 'age'; end
    if nargin < 11 || isempty(batch_source_col), batch_source_col = 'batch'; end
    if nargin < 12 || isempty(age_col), age_col = 'age'; end
    if nargin < 13 || isempty(iqm_variance), iqm_variance = 95; end
    if nargin < 14 || isempty(p_thr), p_thr = 0.05; end
    if nargin < 15 || isempty(max_qcs), max_qcs = Inf; end
    if nargin < 16 || isempty(apply_pca), apply_pca = true; end
    if nargin < 17 || isempty(outfilename), outfilename = 'iqm_harmonised.csv'; end

    if isstring(preserve_covars), preserve_covars = cellstr(preserve_covars); end
    if isstring(adjust_covars), adjust_covars = cellstr(adjust_covars); end
    if isstring(reverse_guard_covars), reverse_guard_covars = cellstr(reverse_guard_covars); end
    if isstring(categorical_covars), categorical_covars = cellstr(categorical_covars); end

    if isstring(batch_col), batch_col = char(batch_col); end
    if isstring(subject_col), subject_col = char(subject_col); end
    if isstring(age_source_col), age_source_col = char(age_source_col); end
    if isstring(batch_source_col), batch_source_col = char(batch_source_col); end
    if isstring(age_col), age_col = char(age_col); end

    preserve_covars = unique(cellstr(preserve_covars), 'stable');
    adjust_covars = unique(cellstr(adjust_covars), 'stable');
    reverse_guard_covars = unique(cellstr(reverse_guard_covars), 'stable');
    categorical_covars = unique(cellstr(categorical_covars), 'stable');

    if ~ismember(batch_col, categorical_covars)
        categorical_covars{end+1} = batch_col;
    end

    % Reverse guards only make sense for numeric preserve covariates
    reverse_guard_covars = setdiff(reverse_guard_covars, categorical_covars, 'stable');
    reverse_guard_covars = intersect(reverse_guard_covars, preserve_covars, 'stable');

    %% ---------------------------------------------------------------------
    %% Logging
    %% ---------------------------------------------------------------------
    fprintf('=== IQM Harmonization (Flexible Version) ===\n');
    fprintf('P-value threshold: %.4f\n', p_thr);
    fprintf('Max QCs per volume: %s\n', num2str(max_qcs));
    fprintf('PCA variance threshold: %.1f%%\n', iqm_variance);
    fprintf('Apply PCA: %d\n', apply_pca);
    fprintf('Preserve covars: %s\n', strjoin(preserve_covars, ', '));
    fprintf('Adjust covars: %s\n', strjoin(adjust_covars, ', '));
    fprintf('Categorical covars: %s\n', strjoin(categorical_covars, ', '));
    fprintf('Reverse-guard covars: %s\n\n', strjoin(reverse_guard_covars, ', '));

    %% ---------------------------------------------------------------------
    %% Validate inputs
    %% ---------------------------------------------------------------------
    if ~iscell(qc_list), error('qc_list must be a cell array'); end
    if ~iscell(idp_list), error('idp_list must be a cell array'); end

    missing_idps = ~ismember(idp_list, data.Properties.VariableNames);
    if any(missing_idps)
        error('Missing IDP columns: %s', strjoin(idp_list(missing_idps), ', '));
    end

    missing_qcs = ~ismember(qc_list, data.Properties.VariableNames);
    if any(missing_qcs)
        error('Missing QC columns: %s', strjoin(qc_list(missing_qcs), ', '));
    end

    if ~ismember(subject_col, data.Properties.VariableNames)
        error('Missing subject column: %s', subject_col);
    end
    if ~ismember(age_source_col, data.Properties.VariableNames) && ~strcmp(age_source_col, age_col)
        error('Missing age source column: %s', age_source_col);
    end
    if ~ismember(batch_source_col, data.Properties.VariableNames) && ~strcmp(batch_source_col, batch_col)
        error('Missing batch source column: %s', batch_source_col);
    end

    %% ---------------------------------------------------------------------
    %% Prepare working table
    %% ---------------------------------------------------------------------
    data_model = data;

    if ~ismember(age_col, data_model.Properties.VariableNames)
        data_model.(age_col) = zscore(data_model.(age_source_col));
    end
    if ~ismember(batch_col, data_model.Properties.VariableNames)
        data_model.(batch_col) = data_model.(batch_source_col);
    end

    data_model.(subject_col) = categorical(data_model.(subject_col));
    if ismember(batch_col, data_model.Properties.VariableNames)
        data_model.(batch_col) = categorical(data_model.(batch_col));
    end
    for i = 1:numel(categorical_covars)
        c = categorical_covars{i};
        if ismember(c, data_model.Properties.VariableNames)
            data_model.(c) = categorical(data_model.(c));
        end
    end

    for i = 1:numel(preserve_covars)
        c = preserve_covars{i};
        if ~ismember(c, data_model.Properties.VariableNames)
            error('Missing preserve covariate column: %s', c);
        end
    end
    for i = 1:numel(adjust_covars)
        c = adjust_covars{i};
        if ~ismember(c, data_model.Properties.VariableNames)
            error('Missing adjust covariate column: %s', c);
        end
    end

    %% ---------------------------------------------------------------------
    %% STEP 1: QC basis
    %% ---------------------------------------------------------------------
    fprintf('Step 1: Preparing QC metrics...\n');

    qc_matrix = table2array(data_model(:, qc_list));

    std_tmp_qc = std(qc_matrix, 0, 1)';
    remove_qc = find(std_tmp_qc == 0);
    if ~isempty(remove_qc)
        fprintf('  Removing %d zero-variance QCs\n', length(remove_qc));
        fprintf('  %s\n', strjoin(qc_list(remove_qc), ', '));
        qc_matrix(:, remove_qc) = [];
        qc_list_clean = qc_list;
        qc_list_clean(remove_qc) = [];
    else
        qc_list_clean = qc_list;
    end

    if isempty(qc_matrix)
        error('No QC variables left after removing zero-variance columns.');
    end

    qc_z = zscore(qc_matrix);

    if apply_pca
        [~, scores_iqm, ~, ~, explained_iqm] = pca(qc_z);
        cumulative_explained_iqm = cumsum(explained_iqm);
        num_qc_features = find(cumulative_explained_iqm >= iqm_variance, 1);

        qc_feature_names = arrayfun(@(i) sprintf('QC%d', i), 1:num_qc_features, 'UniformOutput', false);
        for i = 1:num_qc_features
            data_model.(qc_feature_names{i}) = scores_iqm(:, i);
        end

        qc_basis = 'PCA';
        fprintf('  PCA input shape: (%d, %d)\n', size(qc_z,1), size(qc_z,2));
        fprintf('  Retaining PCs: %d (%.2f%%)\n', num_qc_features, cumulative_explained_iqm(num_qc_features));
        pca_num_components = num_qc_features;
        pca_variance_explained = cumulative_explained_iqm(num_qc_features);
    else
        qc_feature_names = qc_list_clean;
        num_qc_features = numel(qc_feature_names);
        for i = 1:num_qc_features
            data_model.(qc_feature_names{i}) = qc_z(:, i);
        end

        qc_basis = 'RAW';
        fprintf('  Using raw QC variables directly: %d features\n', num_qc_features);
        pca_num_components = [];
        pca_variance_explained = [];
    end

    %% ---------------------------------------------------------------------
    %% Build detail row template
    %% ---------------------------------------------------------------------
    row_template = struct();
    row_template.volume = '';
    row_template.qc = '';
    row_template.stage = '';
    row_template.p_batch = NaN;
    row_template.p_proxy = NaN;
    row_template.passed = false;
    row_template.fit_family = '';
    row_template.fit_method = '';

    for i = 1:numel(preserve_covars)
        row_template.(['p_preserve_' preserve_covars{i}]) = NaN;
    end
    for i = 1:numel(reverse_guard_covars)
        row_template.(['p_reverse_' reverse_guard_covars{i}]) = NaN;
    end

    add_detail_rows = repmat(row_template, 0, 1);
    mult_detail_rows = repmat(row_template, 0, 1);

    good_pcs_add = zeros(num_qc_features, length(idp_list));
    good_pcs_mult = zeros(num_qc_features, length(idp_list));
    all_out_add = cell(length(idp_list), 2);
    all_out_mult = cell(length(idp_list), 2);

    %% ---------------------------------------------------------------------
    %% STEP 2: ADDITIVE QC SELECTION
    %% ---------------------------------------------------------------------
    fprintf('\nStep 2: Selecting QCs for ADDITIVE correction...\n');
    fprintf('  Criteria: batch-driven, not preserve-covariate-driven\n');

    for v = 1:length(idp_list)
        volname = idp_list{v};
        fprintf('\n=== Additive selection: %s ===\n', volname);

        rows = repmat(row_template, num_qc_features, 1);

        for qc_i = 1:num_qc_features
            qcname = qc_feature_names{qc_i};
            row = row_template;
            row.volume = volname;
            row.qc = qcname;
            row.stage = 'additive';

            try
                % one stable model per QC, then ANOVA
                base_terms = [preserve_covars, adjust_covars, {batch_col, volname}];
                base_terms = unique(base_terms, 'stable');
                rhs = strjoin(base_terms, ' + ');
                formula_qc = sprintf('%s ~ %s + (1|%s)', qcname, rhs, subject_col);

                mdl_qc = fitlme(data_model, formula_qc);
                tbl_qc = anova(mdl_qc);

                preserve_pass = true;
                for c = 1:numel(preserve_covars)
                    cov = preserve_covars{c};
                    idx = find(strcmp(tbl_qc.Term, cov), 1);
                    if isempty(idx)
                        p_cov = NaN;
                    else
                        p_cov = tbl_qc.pValue(idx);
                    end
                    row.(['p_preserve_' cov]) = p_cov;
                    preserve_pass = preserve_pass && (p_cov > p_thr);
                end

                idx = find(strcmp(tbl_qc.Term, batch_col), 1);
                if isempty(idx)
                    p_batch = NaN;
                else
                    p_batch = tbl_qc.pValue(idx);
                end
                row.p_batch = p_batch;

                idx = find(strcmp(tbl_qc.Term, volname), 1);
                if isempty(idx)
                    p_proxy = NaN;
                else
                    p_proxy = tbl_qc.pValue(idx);
                end
                row.p_proxy = p_proxy;

                reverse_pass = true;
                for c = 1:numel(reverse_guard_covars)
                    cov = reverse_guard_covars{c};
                    other_preserve = setdiff(preserve_covars, {cov}, 'stable');
                    rev_terms = [{qcname}, other_preserve, adjust_covars, {batch_col, volname}];
                    rev_terms = unique(rev_terms, 'stable');
                    rev_rhs = strjoin(rev_terms, ' + ');
                    rev_formula = sprintf('%s ~ %s + (1|%s)', cov, rev_rhs, subject_col);

                    mdl_rev = fitlme(data_model, rev_formula);
                    tbl_rev = anova(mdl_rev);

                    idx2 = find(strcmp(tbl_rev.Term, qcname), 1);
                    if isempty(idx2)
                        p_rev = NaN;
                    else
                        p_rev = tbl_rev.pValue(idx2);
                    end
                    row.(['p_reverse_' cov]) = p_rev;
                    reverse_pass = reverse_pass && (p_rev > p_thr);
                end

                pass_filter = preserve_pass && (p_batch < p_thr) && (p_proxy > p_thr) && reverse_pass;
                row.passed = pass_filter;
                good_pcs_add(qc_i, v) = pass_filter;

            catch ME
                warning('Error in additive selection for %s, %s: %s', volname, qcname, ME.message);
                row.passed = false;
                good_pcs_add(qc_i, v) = false;
            end

            rows(qc_i) = row;
            add_detail_rows(end+1, 1) = row; %#ok<AGROW>
        end

        all_out_add{v, 1} = volname;
        all_out_add{v, 2} = rows;
    end

    n_qcs_add = sum(good_pcs_add, 1);
    fprintf('  QCs selected: min=%d, max=%d, mean=%.1f\n', ...
        min(n_qcs_add), max(n_qcs_add), mean(n_qcs_add));

    %% ---------------------------------------------------------------------
    %% STEP 3: MULTIPLICATIVE QC SELECTION
    %% ---------------------------------------------------------------------
    fprintf('\nStep 3: Selecting QCs for MULTIPLICATIVE correction...\n');
    fprintf('  KEY: Using ORIGINAL volumes (not additive-corrected) for variance proxy\n');
    fprintf('  Criteria: batch-driven, NOT strongly variance-driven\n');

    for v = 1:length(idp_list)
        volname = idp_list{v};
        fprintf('\n=== Multiplicative selection: %s ===\n', volname);

        rows = repmat(row_template, num_qc_features, 1);

        % biology-only residuals from original volume
        bio_terms = [preserve_covars, adjust_covars];
        bio_terms = unique(bio_terms, 'stable');
        bio_rhs = strjoin(bio_terms, ' + ');
        bio_formula = sprintf('%s ~ %s + (1|%s)', volname, bio_rhs, subject_col);

        lme_bio = fitlme(data_model, bio_formula);
        res = residuals(lme_bio);
        data_model.log_r2_tmp = log(max(res.^2, eps));

        for qc_i = 1:num_qc_features
            qcname = qc_feature_names{qc_i};
            row = row_template;
            row.volume = volname;
            row.qc = qcname;
            row.stage = 'multiplicative';

            try
                base_terms = [preserve_covars, adjust_covars, {batch_col, 'log_r2_tmp'}];
                base_terms = unique(base_terms, 'stable');
                rhs = strjoin(base_terms, ' + ');
                formula_var = sprintf('%s ~ %s + (1|%s)', qcname, rhs, subject_col);

                mdl_var = fitlme(data_model, formula_var);
                tbl_var = anova(mdl_var);

                preserve_pass = true;
                for c = 1:numel(preserve_covars)
                    cov = preserve_covars{c};
                    idx = find(strcmp(tbl_var.Term, cov), 1);
                    if isempty(idx)
                        p_cov = NaN;
                    else
                        p_cov = tbl_var.pValue(idx);
                    end
                    row.(['p_preserve_' cov]) = p_cov;
                    preserve_pass = preserve_pass && (p_cov > p_thr);
                end

                idx = find(strcmp(tbl_var.Term, batch_col), 1);
                if isempty(idx)
                    p_batch = NaN;
                else
                    p_batch = tbl_var.pValue(idx);
                end
                row.p_batch = p_batch;

                idx = find(strcmp(tbl_var.Term, 'log_r2_tmp'), 1);
                if isempty(idx)
                    p_proxy = NaN;
                else
                    p_proxy = tbl_var.pValue(idx);
                end
                row.p_proxy = p_proxy;

                reverse_pass = true;
                for c = 1:numel(reverse_guard_covars)
                    cov = reverse_guard_covars{c};
                    other_preserve = setdiff(preserve_covars, {cov}, 'stable');
                    rev_terms = [{qcname}, other_preserve, adjust_covars, {batch_col, 'log_r2_tmp'}];
                    rev_terms = unique(rev_terms, 'stable');
                    rev_rhs = strjoin(rev_terms, ' + ');
                    rev_formula = sprintf('%s ~ %s + (1|%s)', cov, rev_rhs, subject_col);

                    mdl_rev = fitlme(data_model, rev_formula);
                    tbl_rev = anova(mdl_rev);

                    idx2 = find(strcmp(tbl_rev.Term, qcname), 1);
                    if isempty(idx2)
                        p_rev = NaN;
                    else
                        p_rev = tbl_rev.pValue(idx2);
                    end
                    row.(['p_reverse_' cov]) = p_rev;
                    reverse_pass = reverse_pass && (p_rev > p_thr);
                end

                pass_filter = preserve_pass && (p_batch < p_thr) && (p_proxy > p_thr) && reverse_pass;
                row.passed = pass_filter;
                good_pcs_mult(qc_i, v) = pass_filter;

            catch ME
                warning('Error in multiplicative selection for %s, %s: %s', volname, qcname, ME.message);
                row.passed = false;
                good_pcs_mult(qc_i, v) = false;
            end

            rows(qc_i) = row;
            mult_detail_rows(end+1, 1) = row; %#ok<AGROW>
        end

        all_out_mult{v, 1} = volname;
        all_out_mult{v, 2} = rows;

        if ismember('log_r2_tmp', data_model.Properties.VariableNames)
            data_model.log_r2_tmp = [];
        end
    end

    n_qcs_mult = sum(good_pcs_mult, 1);
    fprintf('  QCs selected: min=%d, max=%d, mean=%.1f\n', ...
        min(n_qcs_mult), max(n_qcs_mult), mean(n_qcs_mult));

    %% ---------------------------------------------------------------------
    %% STEP 4: APPLY CORRECTIONS
    %% ---------------------------------------------------------------------
    fprintf('\nStep 4: Applying corrections...\n');

    qc_selection = struct();
    qc_selection.volumes = idp_list;
    qc_selection.qc_basis = qc_basis;
    qc_selection.apply_pca = apply_pca;
    qc_selection.additive_qcs = cell(length(idp_list), 1);
    qc_selection.multiplicative_qcs = cell(length(idp_list), 1);
    qc_selection.additive_count = zeros(length(idp_list), 1);
    qc_selection.multiplicative_count = zeros(length(idp_list), 1);
    qc_selection.multiplicative_model_pval = nan(length(idp_list), 1);
    qc_selection.multiplicative_applied = false(length(idp_list), 1);
    qc_selection.p_threshold = p_thr;
    qc_selection.max_qcs = max_qcs;
    qc_selection.preserve_covars = preserve_covars;
    qc_selection.adjust_covars = adjust_covars;
    qc_selection.reverse_guard_covars = reverse_guard_covars;
    qc_selection.categorical_covars = categorical_covars;
    qc_selection.pca_info = struct();
    qc_selection.pca_info.num_components = pca_num_components;
    qc_selection.pca_info.variance_explained = pca_variance_explained;
    qc_selection.pca_info.qc_feature_names = qc_feature_names;

    y_harmonised = nan(size(data_model, 1), length(idp_list));

    for v = 1:length(idp_list)
        volname = idp_list{v};
        fprintf('\n=== Volume: %s ===\n', volname);

        %% ADDITIVE CORRECTION
        qc_mask_add = find(good_pcs_add(:, v) == 1);
        if ~isinf(max_qcs) && length(qc_mask_add) > max_qcs
            qc_mask_add = qc_mask_add(1:max_qcs);
        end

        if isempty(qc_mask_add)
            fprintf('  Additive correction: skipped (no QCs passed)\n');
            tmp_noAdditive = data_model.(volname);
            qc_selection.additive_qcs{v} = {};
            qc_selection.additive_count(v) = 0;
        else
            qc_vars_add = qc_feature_names(qc_mask_add);
            fprintf('  Additive correction: applying %d QC(s) -> %s\n', ...
                length(qc_vars_add), strjoin(qc_vars_add, ', '));

            additive_terms = [preserve_covars, adjust_covars, qc_vars_add];
            additive_terms = unique(additive_terms, 'stable');
            rhs = strjoin(additive_terms, ' + ');
            formula_add = sprintf('%s ~ %s + (1|%s)', volname, rhs, subject_col);

            lme_add = fitlme(data_model, formula_add);

            beta_add = fixedEffects(lme_add);
            beta_names_add = lme_add.CoefficientNames;
            idx_qc_add = ismember(beta_names_add, qc_vars_add);
            X_qc_add = data_model{:, beta_names_add(idx_qc_add)};
            additive_effect = X_qc_add * beta_add(idx_qc_add);

            tmp_noAdditive = data_model.(volname) - additive_effect;
            qc_selection.additive_qcs{v} = qc_vars_add;
            qc_selection.additive_count(v) = length(qc_mask_add);
        end

        tmp_noAdditive = max(tmp_noAdditive, 1e-6);

        %% MULTIPLICATIVE CORRECTION
        qc_mask_mult = find(good_pcs_mult(:, v) == 1);
        if ~isinf(max_qcs) && length(qc_mask_mult) > max_qcs
            qc_mask_mult = qc_mask_mult(1:max_qcs);
        end

        if isempty(qc_mask_mult)
            fprintf('  Multiplicative correction: skipped (no QCs passed)\n');
            y_h = tmp_noAdditive;
            qc_selection.multiplicative_qcs{v} = {};
            qc_selection.multiplicative_count(v) = 0;
            qc_selection.multiplicative_applied(v) = false;
            qc_selection.multiplicative_model_pval(v) = NaN;
        else
            qc_vars_mult = qc_feature_names(qc_mask_mult);
            fprintf('  Multiplicative selection: %d QC(s) passed -> %s\n', ...
                length(qc_vars_mult), strjoin(qc_vars_mult, ', '));

            data_model.tmp_noAdditive = tmp_noAdditive;
            data_model.log_tmp_noAdditive = log(data_model.tmp_noAdditive);

            red_terms = [preserve_covars, adjust_covars];
            red_terms = unique(red_terms, 'stable');
            red_rhs = strjoin(red_terms, ' + ');
            formula_red = sprintf('log_tmp_noAdditive ~ %s + (1|%s)', red_rhs, subject_col);

            full_terms = [preserve_covars, adjust_covars, qc_vars_mult];
            full_terms = unique(full_terms, 'stable');
            full_rhs = strjoin(full_terms, ' + ');
            formula_full = sprintf('log_tmp_noAdditive ~ %s + (1|%s)', full_rhs, subject_col);

            try
                lme_red = fitlme(data_model, formula_red);
                lme_full = fitlme(data_model, formula_full);
                cmp = compare(lme_red, lme_full);
                model_pval = cmp.pValue(2);
                qc_selection.multiplicative_model_pval(v) = model_pval;

                fprintf('  Multiplicative model comparison p-value: %.6g\n', model_pval);

                if model_pval < p_thr
                    beta_mult = fixedEffects(lme_full);
                    beta_names_mult = lme_full.CoefficientNames;
                    idx_qc_mult = ismember(beta_names_mult, qc_vars_mult);
                    X_qc_mult = data_model{:, beta_names_mult(idx_qc_mult)};

                    lp = X_qc_mult * beta_mult(idx_qc_mult);
                    lp_centered = lp - mean(lp, 'omitnan');
                    multiplicative_effect = exp(lp_centered);

                    y_h = tmp_noAdditive ./ multiplicative_effect;
                    qc_selection.multiplicative_qcs{v} = qc_vars_mult;
                    qc_selection.multiplicative_count(v) = length(qc_mask_mult);
                    qc_selection.multiplicative_applied(v) = true;
                    fprintf('  Multiplicative correction: applied\n');
                else
                    y_h = tmp_noAdditive;
                    qc_selection.multiplicative_qcs{v} = {};
                    qc_selection.multiplicative_count(v) = 0;
                    qc_selection.multiplicative_applied(v) = false;
                    fprintf('  Multiplicative correction: skipped (model p >= threshold)\n');
                end
            catch ME
                warning('Model comparison failed for %s: %s', volname, ME.message);
                y_h = tmp_noAdditive;
                qc_selection.multiplicative_qcs{v} = {};
                qc_selection.multiplicative_count(v) = 0;
                qc_selection.multiplicative_applied(v) = false;
                qc_selection.multiplicative_model_pval(v) = NaN;
            end

            if ismember('log_tmp_noAdditive', data_model.Properties.VariableNames)
                data_model.log_tmp_noAdditive = [];
            end
            if ismember('tmp_noAdditive', data_model.Properties.VariableNames)
                data_model.tmp_noAdditive = [];
            end
        end

        data.(['harmonised_', volname]) = y_h;
        y_harmonised(:, v) = y_h;
    end

    %% ---------------------------------------------------------------------
    %% STEP 5: SAVE OUTPUTS
    %% ---------------------------------------------------------------------
    fprintf('\nStep 5: Saving results...\n');

    writetable(data, outfilename);
    fprintf('  Data saved: %s\n', outfilename);

    summary_tbl = table( ...
        qc_selection.volumes, ...
        qc_selection.additive_count, ...
        qc_selection.multiplicative_count, ...
        qc_selection.multiplicative_applied, ...
        qc_selection.multiplicative_model_pval, ...
        'VariableNames', {'volume', 'additive_count', 'multiplicative_count', 'multiplicative_applied', 'multiplicative_model_pval'});

    writetable(summary_tbl, 'qc_selection_summary.csv');
    fprintf('  Summary saved: qc_selection_summary.csv\n');

    if isempty(add_detail_rows)
        add_detail_tbl = table();
    else
        add_detail_tbl = struct2table(add_detail_rows);
    end
    writetable(add_detail_tbl, 'qc_selection_additive_details.csv');
    fprintf('  Additive details saved: qc_selection_additive_details.csv\n');

    if isempty(mult_detail_rows)
        mult_detail_tbl = table();
    else
        mult_detail_tbl = struct2table(mult_detail_rows);
    end
    writetable(mult_detail_tbl, 'qc_selection_multiplicative_details.csv');
    fprintf('  Multiplicative details saved: qc_selection_multiplicative_details.csv\n');

    selected_additive = cell(length(qc_selection.additive_qcs), 1);
    for i = 1:numel(qc_selection.additive_qcs)
        q = qc_selection.additive_qcs{i};
        if isempty(q)
            selected_additive{i} = '';
        else
            selected_additive{i} = strjoin(q, ';');
        end
    end

    selected_multiplicative = cell(length(qc_selection.multiplicative_qcs), 1);
    for i = 1:numel(qc_selection.multiplicative_qcs)
        q = qc_selection.multiplicative_qcs{i};
        if isempty(q)
            selected_multiplicative{i} = '';
        else
            selected_multiplicative{i} = strjoin(q, ';');
        end
    end

    selected_additive_tbl = table(string(qc_selection.volumes(:)), string(selected_additive(:)), ...
        'VariableNames', {'volume', 'selected_qcs'});
    writetable(selected_additive_tbl, 'selected_additive_qcs_by_volume.csv');

    selected_multiplicative_tbl = table(string(qc_selection.volumes(:)), string(selected_multiplicative(:)), ...
        'VariableNames', {'volume', 'selected_qcs'});
    writetable(selected_multiplicative_tbl, 'selected_multiplicative_qcs_by_volume.csv');

    fprintf('  Selected QC lists saved.\n');

    save('qc_selection_working.mat', ...
        'qc_selection', 'all_out_add', 'all_out_mult', ...
        'good_pcs_add', 'good_pcs_mult', ...
        'summary_tbl', 'add_detail_tbl', 'mult_detail_tbl');

    fprintf('  QC selection saved: qc_selection_working.mat\n');

    %% Summary printout
    fprintf('\n=== SUMMARY ===\n');
    fprintf('Additive correction:\n');
    fprintf('  Applied to: %d/%d volumes\n', sum(qc_selection.additive_count > 0), length(idp_list));

    if sum(qc_selection.additive_count > 0) > 0
        idx = qc_selection.additive_count > 0;
        fprintf('  QCs used: min=%d, max=%d, mean=%.1f\n', ...
            min(qc_selection.additive_count(idx)), ...
            max(qc_selection.additive_count(idx)), ...
            mean(qc_selection.additive_count(idx)));
    end

    fprintf('\nMultiplicative correction:\n');
    fprintf('  QCs selected for: %d/%d volumes\n', sum(qc_selection.multiplicative_count > 0), length(idp_list));
    fprintf('  Actually applied: %d/%d volumes\n', sum(qc_selection.multiplicative_applied), length(idp_list));

    if sum(qc_selection.multiplicative_applied) > 0
        idx = qc_selection.multiplicative_applied;
        fprintf('  QCs used: min=%d, max=%d, mean=%.1f\n', ...
            min(qc_selection.multiplicative_count(idx)), ...
            max(qc_selection.multiplicative_count(idx)), ...
            mean(qc_selection.multiplicative_count(idx)));
        fprintf('  Model p-values: min=%.4f, max=%.4f, median=%.4f\n', ...
            min(qc_selection.multiplicative_model_pval(idx)), ...
            max(qc_selection.multiplicative_model_pval(idx)), ...
            median(qc_selection.multiplicative_model_pval(idx)));
    end

    fprintf('\n=== KEY FEATURES ===\n');
    fprintf('✓ QCs selected from %s basis\n', qc_basis);
    fprintf('✓ Preserve covariates screened in forward QC models\n');
    fprintf('✓ Reverse-guard applied only to numeric preserve covariates\n');
    fprintf('✓ Additive correction always attempted\n');
    fprintf('✓ Multiplicative correction gated by model comparison\n');
    fprintf('\nDone!\n');
end