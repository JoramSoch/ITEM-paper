% Project directories script
% _
% This script generates project directories for the current study.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% Date  : 10/11/2017, 12:50


% set directories
stud_dir = 'C:\Joram\projects\VisRec\DataSets\Heinzle_et_al_2011\';
bids_dir = strcat(stud_dir,'BIDS\');
tool_dir = strcat(stud_dir,'tools\');

% save directories
save('project_directories.mat','stud_dir','bids_dir','tool_dir');