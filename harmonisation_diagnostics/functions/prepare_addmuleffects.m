function prepare_addmuleffects(yourStruct, idpNames, outputDir)
% createIDPCSVs - Create CSV files for each field with base columns from Raw and IDPs
%
% Syntax: createIDPCSVs(yourStruct, idpNames, outputDir)
%
% Inputs:
%   yourStruct - Structure with fields containing tables
%   idpNames   - Cell array of IDP column names to extract
%   outputDir  - (Optional) Directory to save CSV files. Default is current directory
%
% Example:
%   idpNames = {'Volume', 'Thickness', 'Area'};
%   createIDPCSVs(myStruct, idpNames, './output/');

% Handle optional output directory
if nargin < 3
    outputDir = './';
end

% Ensure output directory ends with separator
if ~endsWith(outputDir, filesep)
    outputDir = [outputDir filesep];
end

% Create output directory if it doesn't exist
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% Configuration
baseCols = {'subjectID', 'timepoint', 'Site', 'Scanner', 'Final_Age'};
fields = fieldnames(yourStruct);

% First, check if subjectID and timepoint match across all fields
fprintf('=== Checking subjectID and timepoint consistency ===\n');
subjectIDs = cellfun(@(f) yourStruct.(f).subjectID, fields, 'UniformOutput', false);
timepoints = cellfun(@(f) yourStruct.(f).timepoint, fields, 'UniformOutput', false);

subjectIDMatch = all(cellfun(@(x) isequal(x, subjectIDs{1}), subjectIDs));
timepointMatch = all(cellfun(@(x) isequal(x, timepoints{1}), timepoints));

if subjectIDMatch && timepointMatch
    fprintf('✓ All subjectID and timepoint columns match!\n\n');
else
    if ~subjectIDMatch
        error('✗ subjectID columns do NOT match across fields!');
    end
    if ~timepointMatch
        error('✗ timepoint columns do NOT match across fields!');
    end
end

% Get base columns from Raw
fprintf('=== Creating CSVs ===\n');
baseTable = yourStruct.Raw(:, baseCols);
baseTable.Properties.VariableNames{'Final_Age'} = 'age';
baseTable.zscore_age = zscore(baseTable.age);

% Process each field
for i = 1:length(fields)
    fieldName = fields{i};
    currentTable = yourStruct.(fieldName);
    
    % Start with base columns from Raw
    outputTable = baseTable;
    
    % Add each IDP from current field
    for j = 1:length(idpNames)
        idp = idpNames{j};
        colNames = currentTable.Properties.VariableNames;
        
        % Search for harmonised version first
        harmonisedPattern = ['^harmonised.*' idp];
        harmonisedMatch = ~cellfun(@isempty, regexpi(colNames, harmonisedPattern));
        
        if any(harmonisedMatch)
            % Harmonised version found - keep the harmonised name
            matchedCol = colNames{find(harmonisedMatch, 1)};
            outputTable.(matchedCol) = currentTable.(matchedCol);
            fprintf('  Field "%s", IDP "%s": Using "%s"\n', fieldName, idp, matchedCol);
        elseif ismember(idp, colNames)
            % Original IDP found
            % For non-Raw fields, rename column to harmonised_<idp>
            if ~strcmp(fieldName, 'Raw')
                outputTable.(['harmonised_' idp]) = currentTable.(idp);
                fprintf('  Field "%s", IDP "%s": Using original (renamed to harmonised_%s)\n', fieldName, idp, idp);
            else
                outputTable.(idp) = currentTable.(idp);
                fprintf('  Field "%s", IDP "%s": Using original\n', fieldName, idp);
            end
        else
            fprintf('  Field "%s", IDP "%s": NOT FOUND - skipping\n', fieldName, idp);
        end
    end
    
    % Write to CSV
    csvFileName = [outputDir fieldName '.csv'];
    writetable(outputTable, csvFileName);
    fprintf('✓ Created: %s\n\n', csvFileName);
end

fprintf('=== All CSVs created successfully ===\n');

end