function harmonised_data = harmonise_LR_LME(data_table, batch_variable, harmonisation_method, idp_names, varargin)
% HARMONISE_BRAIN_VOLUMES Harmonise brain imaging measures by removing batch effects
%
% Inputs:
%   data_table           - Table containing imaging data and covariates
%   batch_variable       - Batch variable(s): 'site', 'scanner', or {'scanner','site'}
%                         Note: These refer to the actual column names specified in
%                         'site_var' and 'scanner_var' parameters
%   harmonisation_method - 'LME' (mixed effects) or 'LR' (linear regression)
%   idp_names           - Cell array of imaging measure column names to harmonise
%   varargin            - Optional name-value pairs:
%                         'age_var' (default: 'zscore_age') - name of age column
%                         'timepoint_var' (default: 'timepoint') - name of timepoint column
%                         'subject_var' (default: 'subjectID') - name of subject ID column
%                         'site_var' (default: 'site') - name of site column in table
%                         'scanner_var' (default: 'scanner') - name of scanner column in table
%                         'zscore_age' (default: true) - whether to z-score age automatically
%
% Output:
%   harmonised_data     - Table with harmonised imaging measures (same size as input)
%
% Examples:
%   % Site as fixed effect with LME
%   harm_data = harmonise_brain_volumes(data, 'site', 'LME', {'volume_1', 'volume_2'});
%
%   % Scanner as fixed effect with LR, custom column names
%   harm_data = harmonise_brain_volumes(data, 'scanner', 'LR', {'volume_1'}, ...
%       'scanner_var', 'scanner_type', 'age_var', 'age_years');
%
%   % Scanner fixed + site random with LME
%   harm_data = harmonise_brain_volumes(data, {'scanner','site'}, 'LME', {'volume_1'}, ...
%       'site_var', 'acquisition_site', 'scanner_var', 'mri_scanner');

    % Parse optional inputs
    p = inputParser;
    addParameter(p, 'age_var', 'zscore_age', @(x) ischar(x) || isstring(x));
    addParameter(p, 'timepoint_var', 'timepoint', @(x) ischar(x) || isstring(x));
    addParameter(p, 'subject_var', 'subjectID', @(x) ischar(x) || isstring(x));
    addParameter(p, 'site_var', 'site', @(x) ischar(x) || isstring(x));
    addParameter(p, 'scanner_var', 'scanner', @(x) ischar(x) || isstring(x));
    addParameter(p, 'zscore_age', true, @(x) islogical(x) && isscalar(x));
    parse(p, varargin{:});
    
    age_var = p.Results.age_var;
    timepoint_var = p.Results.timepoint_var;
    subject_var = p.Results.subject_var;
    site_var = p.Results.site_var;
    scanner_var = p.Results.scanner_var;
    do_zscore_age = p.Results.zscore_age;
    
    % Z-score age if requested
    if do_zscore_age && ismember(age_var, data_table.Properties.VariableNames)
        if isnumeric(data_table.(age_var))
            data_table.(age_var) = zscore(data_table.(age_var));
            fprintf('Z-scored age variable: %s\n', age_var);
        end
    end
    
    % Validate inputs
    validateInputs(data_table, batch_variable, harmonisation_method, idp_names, ...
                   age_var, timepoint_var, subject_var, site_var, scanner_var);
    
    % Initialize output table
    harmonised_data = data_table;
    
    % Determine batch configuration
    if ischar(batch_variable)
        batch_variable = {batch_variable};
    end
    
    % Map batch_variable to actual column names
    actual_batch_vars = cell(size(batch_variable));
    for i = 1:length(batch_variable)
        if strcmp(batch_variable{i}, 'site')
            actual_batch_vars{i} = site_var;
        elseif strcmp(batch_variable{i}, 'scanner')
            actual_batch_vars{i} = scanner_var;
        else
            error('batch_variable must contain ''site'' and/or ''scanner''');
        end
    end
    
    % Check if site is included as random effect
    has_site_random = length(batch_variable) == 2 && ...
                      any(strcmp(batch_variable, 'site')) && ...
                      any(strcmp(batch_variable, 'scanner'));
    
    % Build formula components
    if has_site_random
        % Scanner fixed, site random
        fixed_batch = {scanner_var};
        random_batch = {site_var};
    else
        % All batch variables as fixed effects
        fixed_batch = actual_batch_vars;
        random_batch = {};
    end
    
    % Harmonise each imaging measure
    fprintf('Harmonising %d imaging measures using %s method...\n', ...
            length(idp_names), harmonisation_method);
    
    for ii = 1:length(idp_names)
        idp_name = idp_names{ii};
        
        if mod(ii, 10) == 0
            fprintf('  Processing %d/%d: %s\n', ii, length(idp_names), idp_name);
        end
        
        try
            if strcmp(harmonisation_method, 'LME')
                harmonised_data.(idp_name) = harmonise_lme(data_table, idp_name, ...
                    fixed_batch, random_batch, age_var, timepoint_var, subject_var);
            elseif strcmp(harmonisation_method, 'LR')
                harmonised_data.(idp_name) = harmonise_lr(data_table, idp_name, ...
                    fixed_batch, age_var, timepoint_var);
            end
        catch ME
            warning('Failed to harmonise %s: %s', idp_name, ME.message);
            % Keep original data if harmonisation fails
            harmonised_data.(idp_name) = data_table.(idp_name);
        end
    end
    
    fprintf('Harmonisation complete!\n');
end


function harmonised = harmonise_lme(data_table, idp_name, fixed_batch, random_batch, ...
                                     age_var, timepoint_var, subject_var)
% Harmonise using Linear Mixed Effects model
    
    % Build formula
    formula_str = sprintf('%s ~ %s + %s', idp_name, age_var, timepoint_var);
    
    % Add fixed batch effects
    for bb = 1:length(fixed_batch)
        formula_str = [formula_str, ' + ', fixed_batch{bb}];
    end
    
    % Add random effects
    formula_str = [formula_str, sprintf(' + (1|%s)', subject_var)];
    
    if ~isempty(random_batch)
        for bb = 1:length(random_batch)
            formula_str = [formula_str, sprintf(' + (1|%s)', random_batch{bb})];
        end
    end
    
    % Fit model
    mdl = fitlme(data_table, formula_str);
    
    % Get fixed effects design matrix
    X_full = designMatrix(mdl, 'Fixed');
    fixedEff = fixedEffects(mdl);
    varNames = mdl.CoefficientNames;
    
    % Identify batch fixed effects (exclude intercept, age, timepoint)
    batchIdx = ~(contains(varNames, {'(Intercept)', age_var, timepoint_var}));
    
    % Remove fixed batch effects
    batch_contribution = 0;
    if any(batchIdx)
        X_batch = X_full(:, batchIdx);
        batchEff = fixedEff(batchIdx);
        batch_contribution = X_batch * batchEff;
    end
    
    % Get random effects for batch variables (if any)
    random_contribution = 0;
    if ~isempty(random_batch)
        [~, ~, stats] = randomEffects(mdl);
        
        % Extract random effects table
        re_table = stats;
        
        for bb = 1:length(random_batch)
            batch_name = random_batch{bb};
            
            % Find rows corresponding to this random effect
            re_group_idx = strcmp(re_table.Group, batch_name);
            
            if any(re_group_idx)
                re_levels = re_table.Level(re_group_idx);
                re_estimates = re_table.Estimate(re_group_idx);
                
                % Map random effects to observations
                batch_levels = data_table.(batch_name);
                
                % Create mapping
                re_map = containers.Map(re_levels, re_estimates);
                
                % Get random effect for each observation
                batch_re = zeros(height(data_table), 1);
                for jj = 1:height(data_table)
                    level_key = batch_levels(jj);
                    if iscell(level_key)
                        level_key = level_key{1};
                    end
                    level_str = char(string(level_key));
                    
                    if isKey(re_map, level_str)
                        batch_re(jj) = re_map(level_str);
                    end
                end
                
                random_contribution = random_contribution + batch_re;
            end
        end
    end
    
    % Harmonised data: original - fixed batch effects - random batch effects
    harmonised = data_table.(idp_name) - batch_contribution - random_contribution;
end


function harmonised = harmonise_lr(data_table, idp_name, fixed_batch, age_var, timepoint_var)
% Harmonise using Linear Regression
    
    % Build design matrix
    y = data_table.(idp_name);
    
    % Start with age and timepoint
    predictors = {age_var, timepoint_var};
    
    % Add batch variables
    predictors = [predictors, fixed_batch];
    
    % Create design matrix
    X = [];
    predictor_names = {};
    
    for pp = 1:length(predictors)
        var_name = predictors{pp};
        var_data = data_table.(var_name);
        
        if isnumeric(var_data)
            % Numeric variable - add as is
            X = [X, var_data];
            predictor_names{end+1} = var_name;
        else
            % Categorical variable - create dummy variables
            if iscell(var_data) || isstring(var_data)
                var_data = categorical(var_data);
            end
            
            dummy = dummyvar(var_data);
            % Remove first column (reference category)
            dummy = dummy(:, 2:end);
            X = [X, dummy];
            
            % Add names for each dummy variable
            cats = categories(var_data);
            for cc = 2:length(cats)
                predictor_names{end+1} = [var_name, '_', char(cats{cc})];
            end
        end
    end
    
    % Add intercept
    X = [ones(size(X, 1), 1), X];
    predictor_names = [{'Intercept'}, predictor_names];
    
    % Fit linear regression
    beta = X \ y;
    
    % Identify batch effects (exclude intercept, age, timepoint)
    batch_idx = ~(contains(predictor_names, {'Intercept', age_var, timepoint_var}));
    
    % Remove batch effects
    X_batch = X(:, batch_idx);
    beta_batch = beta(batch_idx);
    batch_contribution = X_batch * beta_batch;
    
    harmonised = y - batch_contribution;
end


function validateInputs(data_table, batch_variable, harmonisation_method, idp_names, ...
                       age_var, timepoint_var, subject_var, site_var, scanner_var)
% Validate all inputs
    
    % Check data_table
    assert(istable(data_table), 'data_table must be a MATLAB table');
    
    % Check batch_variable
    if ischar(batch_variable)
        batch_variable = {batch_variable};
    end
    
    valid_batch = {'site', 'scanner'};
    for bb = 1:length(batch_variable)
        assert(ismember(batch_variable{bb}, valid_batch), ...
               'batch_variable must be ''site'', ''scanner'', or {''scanner'',''site''}');
    end
    
    % Check that actual batch columns exist in data_table
    if any(strcmp(batch_variable, 'site'))
        assert(ismember(site_var, data_table.Properties.VariableNames), ...
               'Site column ''%s'' not found in data_table', site_var);
    end
    if any(strcmp(batch_variable, 'scanner'))
        assert(ismember(scanner_var, data_table.Properties.VariableNames), ...
               'Scanner column ''%s'' not found in data_table', scanner_var);
    end
    
    % Check harmonisation_method
    valid_methods = {'LME', 'LR'};
    assert(ismember(harmonisation_method, valid_methods), ...
           'harmonisation_method must be ''LME'' or ''LR''');
    
    % Check idp_names
    assert(iscell(idp_names), 'idp_names must be a cell array');
    for ii = 1:length(idp_names)
        assert(ismember(idp_names{ii}, data_table.Properties.VariableNames), ...
               'Imaging measure ''%s'' not found in data_table', idp_names{ii});
    end
    
    % Check required covariates
    assert(ismember(age_var, data_table.Properties.VariableNames), ...
           'Age variable ''%s'' not found in data_table', age_var);
    assert(ismember(timepoint_var, data_table.Properties.VariableNames), ...
           'Timepoint variable ''%s'' not found in data_table', timepoint_var);
    
    % Check subject variable for LME
    if strcmp(harmonisation_method, 'LME')
        assert(ismember(subject_var, data_table.Properties.VariableNames), ...
               'Subject variable ''%s'' not found in data_table (required for LME)', subject_var);
    end
end