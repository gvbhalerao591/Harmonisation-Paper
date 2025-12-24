function MD = get_MD_numerically_stable(T, yourVar)

    features = table2array(T);
    
    % Step 2: Get the site labels and their groups
    [siteGroups, ~, siteIdx] = unique(yourVar, 'stable');  % siteIdx is an N x 1 vector indicating site group
    
    % Step 3: Determine the number of samples per site (assuming equal counts is not required)
    numSites = numel(siteGroups);
    numFeatures = size(features, 2);
    maxSamples = max(histcounts(siteIdx, 1:numSites+1));
    
    % Step 4: Preallocate a 3D matrix, padding with NaNs if sites have unequal sample counts
    X = nan(maxSamples, numFeatures, numSites);
    
    % Step 5: Fill the 3D matrix
    for s = 1:numSites
        idx = find(siteIdx == s);
        n = numel(idx);
        X(1:n, :, s) = features(idx, :);
    end
    MD = calc_MD_numerically_stable(X);
end
