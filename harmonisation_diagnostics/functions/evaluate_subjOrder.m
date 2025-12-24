function results = evaluate_subjOrder(data1, data2, idp_names, nPerm)
% data1     - Table: IDPs from scanner 1 (pre-harmonisation or raw)
% data2     - Table: IDPs from scanner 2 (paired, same subjects)
% idp_names - Cell array of IDP variable names (e.g., {'GM_vol', 'hippo_vol'})
% nPerm     - Number of permutations for null (e.g., 10000)

    if nargin < 4
        nPerm = 10000;
    end
    
    nIDPs = length(idp_names);
    nSubj = height(data1);
    
    results = table('Size', [nIDPs 4], ...
        'VariableTypes', {'string', 'double', 'double', 'double'}, ...
        'VariableNames', {'IDP', 'SpearmanRho', 'NullMeanRho', 'pValue'});
    
    % Generate null distribution (same for all IDPs)
    null_rhos = zeros(nPerm,1);
    rng('default')
    for i = 1:nPerm
        r1 = randperm(nSubj);
        r2 = randperm(nSubj);
        null_rhos(i) = corr(r1', r2', 'Type', 'Spearman');
    end
    
    % For each IDP
    for i = 1:nIDPs
        idp = idp_names{i};
        x = data1.(idp);
        y = data2.(idp);
        real_rho = corr(x, y, 'Type', 'Spearman');
        
        % p-value (proportion of null >= observed)
        pval = mean(abs(null_rhos) >= abs(real_rho));
    
        results.IDP(i)          = idp;
        results.SpearmanRho(i)  = real_rho;
        results.NullMeanRho(i)  = mean(null_rhos);
        results.pValue(i)       = pval;
    end

end
