maindir     = '/Users/psyc1586_admin/GVB_data/harmonisation_work/00_aa_skullImprovement_Oct2025/00_running_for_paper';
remove_subj = {'UCL006'};

%
idp_list = {'T1_SIENAX_periphGM_norm_vol',...
'T1_SIENAX_CSF_norm_vol',...
'T1_SIENAX_GM_norm_vol',...
'T1_SIENAX_WM_norm_vol',...
'T1_SIENAX_brain_norm_vol',...
'T1_FIRST_left_hippocampus',...
'T1_FIRST_right_hippocampus'}';

%% original data
ori_fname   = 'data_ori.mat';
data        = load(fullfile(maindir, ori_fname));
tabname     = char(fieldnames(data));
ori_data    = data.(tabname);
ori_data    = ori_data(~ismember(ori_data.subjectID, remove_subj),:);

%% linear regression
tmp_dir             = '/Users/psyc1586_admin/GVB_data/harmonisation_work/00_aa_skullImprovement_Oct2025/00_running_for_paper/LR_LME_approaches';
lr_site_adjusted    = readtable(fullfile(tmp_dir,"harmonised_LR_site.csv"));
lr_scanner_adjusted = readtable(fullfile(tmp_dir,"harmonised_LR_scanner.csv"));

%% linear mixed effects
lme_site_adjusted        = readtable(fullfile(tmp_dir, "harmonised_LME_site.csv"));
lme_scanner_adjusted     = readtable(fullfile(tmp_dir, "harmonised_LME_scanner.csv"));
lme_scannersite_adjusted = readtable(fullfile(tmp_dir, "harmonised_LME_scannersite.csv"));

%% IQM approach: 
% Site
tmp_dir           = '/Users/psyc1586_admin/GVB_data/harmonisation_work/00_aa_skullImprovement_Oct2025/00_running_for_paper/IQM_approach';
lme_iqms_site_adjusted = readtable(fullfile(tmp_dir, "iqm_harmonised.csv"));

% Scanner
tmp_dir           = '/Users/psyc1586_admin/GVB_data/harmonisation_work/00_aa_skullImprovement_Oct2025/00_running_for_paper/Exploring_IQM_approaches/IQM_scannerOnly';
lme_iqms_scanner_adjusted = readtable(fullfile(tmp_dir, "iqm_harmonised_scannerOnly.csv"));

%% longitudinal combat
tmp_dir                 = '/Users/psyc1586_admin/GVB_data/harmonisation_work/00_aa_skullImprovement_Oct2025/00_running_for_paper/longComBat';
lc_site_adjusted        = readtable(fullfile(tmp_dir,"harmonised_data_Site_20251014_153420.csv"));
lc_scanner_adjusted     = readtable(fullfile(tmp_dir,"harmonised_data_Scanner_20251014_153420.csv"));
lc_scannersite_adjusted = readtable(fullfile(tmp_dir,"harmonised_data_ScannerF_SiteR_20251014_153420"));

%% Histogram matching
hist_match          = load(fullfile(maindir, 'data_histMatched.mat'));
hist_match          = hist_match.final_hist;
hist_match_adjusted = hist_match(~ismember(hist_match.subjectID,remove_subj),:);

%% SynthSR 
sr_method           = load(fullfile(maindir, "data_synthSRharmonised.mat"));
sr_method           = sr_method.final_sr;
sr_method_adjusted  = sr_method(~ismember(sr_method.subjectID,remove_subj),:);

%% ComBat
combat_tool_path = '/Users/psyc1586_admin/Library/CloudStorage/OneDrive-Nexus365/00_NDCN_work/UCL_Pawel_data/tools/ComBatHarmonization-master/Matlab';
addpath(combat_tool_path)

% 1) Site
data_ori               = ori_data;
zscore_age             = zscore(data_ori.Final_Age);
timepoint              = dummyvar(categorical(data_ori.timepoint));
mod                    = [zscore_age, timepoint(:,2)]; % drop one column of categorical covariate: refer documentation GitHub
batch                  = double(categorical(data_ori.Site));
idps                   = table2array(data_ori(:,idp_list))';

combat_adjusted_para   = combat(idps, batch, mod, 1); % setting parameteric
combat_adjusted_para1  = combat_adjusted_para';
combat_adjusted1       = cell2table([data_ori.subjectID, data_ori.Site, ...
                                data_ori.Group, num2cell(data_ori.Final_Age), data_ori.timepoint,...
                                ],'VariableNames',{'subjectID','Site','groupname','age','timepoint'});
combat_adjusted2       = array2table(combat_adjusted_para1, 'VariableNames',strcat('harmonisedCombat_',idp_list)');
combat_site_adjusted   = [combat_adjusted1, data_ori(:,idp_list), combat_adjusted2];

% 2) Scanner
clear batch
batch                          = double(categorical(data_ori.Scanner));
combat_adjusted_para_scanner   = combat(idps, batch, mod, 1);
combat_adjusted_para1_scanner  = combat_adjusted_para_scanner';
combat_adjusted1_scanner       = cell2table([data_ori.subjectID, data_ori.Scanner, data_ori.Site, ...
                                    data_ori.Group, num2cell(data_ori.Final_Age), data_ori.timepoint,...
                                    ],'VariableNames',{'subjectID','Scanner','Site','groupname','age','timepoint'});
combat_adjusted2_scanner        = array2table(combat_adjusted_para1_scanner, 'VariableNames',strcat('harmonisedCombat_',idp_list)');
combat_scanner_adjusted         = [combat_adjusted1_scanner, data_ori(:,idp_list), combat_adjusted2_scanner];

%% 6) CovBat method
tmp_dir = '/Users/psyc1586_admin/GVB_data/harmonisation_work/00_aa_skullImprovement_Oct2025/00_running_for_paper/CovBat_approach';

% 1) Scanner
covbat_scanner_adjusted = readtable(fullfile(tmp_dir, 'covbat_harmonised_Scanner_output.csv'));

% 2) Site
covbat_site_adjusted = readtable(fullfile(tmp_dir, 'covbat_harmonised_Site_output.csv'));

%% Save all
outdir = fullfile(pwd, 'data');
if ~exist(outdir, 'dir')
    mkdir(outdir)
end

data                           = struct();
data.Raw                       = ori_data;
data.LR_Site                   = lr_site_adjusted;
data.LR_Scanner                = lr_scanner_adjusted;
data.ComBat_Site               = combat_site_adjusted;
data.ComBat_Scanner            = combat_scanner_adjusted;
data.CovBat_Site               = covbat_site_adjusted;
data.CovBat_Scanner            = covbat_scanner_adjusted;
data.LME_Site                  = lme_site_adjusted;
data.LME_Scanner               = lme_scanner_adjusted;
data.LME_ScannerFSiteR         = lme_scannersite_adjusted;
data.LME_IQM_Site              = lme_iqms_site_adjusted;
data.LME_IQM_Scanner           = lme_iqms_scanner_adjusted;
data.LongComBat_Site           = lc_site_adjusted;
data.LongComBat_Scanner        = lc_scanner_adjusted;
data.LongComBat_ScannerFSiteR  = lc_scannersite_adjusted;
data.HistMatch                 = hist_match_adjusted;
data.SynthSR                   = sr_method_adjusted;

save(fullfile(outdir,'data.mat'),'data')

  