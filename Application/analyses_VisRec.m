% Step 0: launch SPM
spm fmri
spm_jobman('initcfg');

% Step 1: pre-processing
tool_dir = pwd;
preproc_meta_batch('subjects_all.mat', 'prepare')
preproc_meta_batch('subjects_all.mat', 'perform')

% Step 2: ITEM analysis
cd(strcat(tool_dir,'/VisRec_ITEM/'));
VisRec_ITEM
cd(strcat(tool_dir,'/VisRec_ITEM/'));
VisRec_ITEM_LS_A
cd(strcat(tool_dir,'/VisRec_ITEM/'));
VisRec_ITEM_LS_S