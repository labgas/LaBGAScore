%% control_use_before_def.m
%
% POSITIVE CONTROL for clean/use_before_def.py. Not a runnable analysis.
%
% This file reproduces, in miniature, the failure that motivated the checker:
% an option is READ eight lines above the guard that defines it. In MATLAB this
% dies at the read - "Unrecognized function or variable 'contrast_objects_tag'"
% - and it dies LATE, because in the real script the read sat near the end,
% after the expensive work and before anything was saved.
%
% checkcode passes this file: every line is syntactically perfect.
%
% EXPECTED: use_before_def.py reports contrast_objects_tag,
%           def line 27 AFTER use line 22.
%
% If a change to the checker makes this file pass, the checker is broken.
%
% -------------------------------------------------------------------------
%
% part of the LaBGAScore script-checker positive controls

printhdr(sprintf('saving contrasts%s', contrast_objects_tag));   % <-- USE (line 22)

results_suffix = 'vox_gm';
myscaling_glm  = 'scaled';

if ~exist('contrast_objects_tag','var')                          % <-- DEF (line 27), too late
    contrast_objects_tag = '';
end

savefilenamedata = fullfile(resultsdir, ...
    ['contrast_data_objects' contrast_objects_tag '.mat']);
