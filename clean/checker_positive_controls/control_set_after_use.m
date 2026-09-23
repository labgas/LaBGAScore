%% control_set_after_use.m
%
% POSITIVE CONTROL for clean/set_after_use.py. Not a runnable analysis.
%
% This file reproduces the failure that motivated the checker: the option is
% guarded near the top, CONSUMED a few lines later, and only THEN set to the
% value the author intended. Nothing errors. The output directory is built from
% the default, and the study's own setting does nothing at all.
%
% That is the results_tag case: eleven decoding runs wrote into the same
% untagged folder and overwrote each other's maps. The statistics were fine -
% results are always recomputed - but every saved artefact was lost.
%
% checkcode passes this file, and so does use_before_def.py: the variable IS
% defined before it is used. Only the ordering is wrong.
%
% EXPECTED: set_after_use.py reports results_tag,
%           def line 26, USED line 29, SET line 36.
%
% If a change to the checker makes this file pass, the checker is broken.
%
% -------------------------------------------------------------------------
%
% part of the LaBGAScore script-checker positive controls

if ~exist('results_tag','var'), results_tag = ''; end            % <-- DEF (line 26)

% OUTPUT DIRECTORIES
tdt_resultsdir = fullfile(basedir, 'TDT', results_tag);          % <-- USE (line 29)
if ~exist(tdt_resultsdir, 'dir')
    mkdir(tdt_resultsdir);
end

% ... ~290 lines of other options in the real script ...

results_tag = 'nocombat_all';                                    % <-- SET (line 36), too late
