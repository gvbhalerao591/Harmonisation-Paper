function [sig_sites_n, outs] = pairwise_siteTests(Mdl, groupVar, dataTable, p_thr)
    fixedEffectsNames = Mdl.CoefficientNames;
    siteTerms = fixedEffectsNames(startsWith(fixedEffectsNames, groupVar));
    
    % Extract full list of site categories
    if ~iscategorical(dataTable.(groupVar))
        dataTable.(groupVar) = categorical(dataTable.(groupVar));
    end

    groupCats = categories(dataTable.(groupVar));
    
    % Determine reference category
    dummySiteNames = erase(siteTerms, [groupVar '_']);
    refSite = setdiff(groupCats, dummySiteNames);
    
    fprintf('Reference site: %s\n', string(refSite));
    
    % Combine all sites including reference
    allSites = [refSite, dummySiteNames];
    nSites = numel(allSites);
    
    ii = 1;
    outs = cell(nchoosek(nSites, 2), 4);
    
    fprintf('Pairwise site comparisons (p-values):\n');
    for i = 1:nSites
        for j = i+1:nSites
            C = zeros(1, length(fixedEffectsNames));
            
            % Only include contrast if not reference
            if ~strcmp(allSites{i}, refSite)
                C(strcmp(fixedEffectsNames, [groupVar '_' allSites{i}])) = 1;
            end
            if ~strcmp(allSites{j}, refSite)
                C(strcmp(fixedEffectsNames, [groupVar '_' allSites{j}])) = -1;
            end
            
            p = coefTest(Mdl, C);
            fprintf('%s vs %s: p = %.5f\n', allSites{i}, allSites{j}, p);

            outs{ii,1} = allSites{i};
            outs{ii,2} = allSites{j};
            outs{ii,3} = p;
            if p < p_thr
               tmp(ii,1)  = 1;
               outs{ii,4} = 1;
            else
                tmp(ii,1)  = 0;
                outs{ii,4} = 0;
            end
            ii = ii + 1;
        end
    end
    sig_sites_n = sum(tmp);
end
