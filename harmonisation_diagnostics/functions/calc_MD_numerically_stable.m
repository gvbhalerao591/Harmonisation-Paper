function [MD, num_retained] = calc_MD_numerically_stable(matrix)
    numSites = size(matrix,3);
    numFeatures = size(matrix,2);
    
    allMeans = zeros(numFeatures, numSites);
    tmpCov = zeros(numFeatures, numFeatures);
    
    for sites = 1:numSites
        allMeans(:, sites) = mean(matrix(:,:,sites),'omitnan');
        tmpCov = tmpCov + cov(matrix(:,:,sites),"omitrows");
    end
    
    overallCov = tmpCov./numSites;
    overallMean = mean(allMeans,2,'omitnan');
    
    % Check condition number and use appropriate solver
    condition_num = cond(overallCov);

    if condition_num > 1e15
        fprintf('Using SVD-based pseudoinverse (condition = %.2e)\n', condition_num);
        
        % Use SVD for maximum numerical stability
        [U, S, V] = svd(overallCov);
        s = diag(S);
        
        % More conservative truncation for very ill-conditioned matrices
        tolerance = max(s) * max(size(overallCov)) * eps;
        s_inv = zeros(size(s));
        s_inv(s > tolerance) = 1./s(s > tolerance);
        
        overallCov_pinv = V * diag(s_inv) * U';
        
        MD = zeros(numSites,1);
        for sites = 1:numSites
            term1 = (allMeans(:,sites) - overallMean);
            delta = term1' * overallCov_pinv * term1;
            MD(sites,1) = sqrt(max(delta, 0)); % Ensure non-negative
        end
        
        num_retained = sum(s > tolerance);
        fprintf('Retaining %d of %d singular values\n', num_retained, numel(s));
    
        
    else
        num_retained = 0;
        fprintf('Retaining %d of %d singular values\n', num_retained, numFeatures);
        fprintf('Using standard solver (condition = %.2e)\n', condition_num);
        
        % Standard calculation for well-conditioned matrices
        MD = zeros(numSites,1);
        for sites = 1:numSites
            term1 = (allMeans(:,sites) - overallMean);
            delta = term1' * (overallCov \ term1);
            MD(sites,1) = sqrt(max(delta, 0));
        end
    end
end