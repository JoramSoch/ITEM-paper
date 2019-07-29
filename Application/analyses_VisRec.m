% Step 1: pre-processing
tool_dir = pwd;
preproc_meta_batch('subjects_all.mat', 'prepare')
preproc_meta_batch('subjects_all.mat', 'perform')

% Step 2: ITEM analysis
cd(strcat(pwd,'/VisRec_ITEM/'));
VisRec_ITEM
VisRec_ITEM_LS_A
VisRec_ITEM_LS_S