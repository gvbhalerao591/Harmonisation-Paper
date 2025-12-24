%%
maindir = '/Users/psyc1586_admin/GVB_data/harmonisation_work/00_aa_skullImprovement_Oct2025/run_harmonisation_pipelines';
fname   = 'for_iqm_harmonisation.mat';

load(fullfile(maindir, fname));
load(fullfile(maindir, 'IQM_approach', 'idp_list'));
load(fullfile(maindir, 'IQM_approach', 'iqm_list'));

data            = final_iqm;
data.age        = zscore(data.Final_Age);
data.zscore_age = zscore(data.Final_Age);
data.batch      = data.Site; % change this when using other batch variable e.g., "Scanner"
data.timepoint  = data.timepoint;
data.subjectID  = data.subjectID;

data           = data(~ismember(data.subjectID, {'UCL006'}),:);
harmonisedData = iqm_harmonise(data, idp_list, iqm_list, 95, 0.05,Inf,'iqm_harmonised.csv');


