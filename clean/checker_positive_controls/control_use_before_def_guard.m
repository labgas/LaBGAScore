%% control_use_before_def_guard.m
%
% SECOND POSITIVE CONTROL for clean/use_before_def.py. Not a runnable analysis.
%
% The sibling control, control_use_before_def.m, spells its default over three
% lines, so the assignment sits on a line of its own and the checker's
% bare-assignment pattern sees it. THIS file uses the ONE-LINE guard idiom
%
%     if ~exist('x','var'), x = []; end
%
% which is the dominant form in the real templates - in prep_3a all eleven
% guarded options are written this way. Until 2026-10-05 the checker built its
% set of option-ish names from bare assignments only, so every option in that
% form recorded ZERO uses and passed unconditionally: on prep_3a the check could
% only ever return PASS, and it did, while cv_seed_mvpa_reg_cov was being read
% four lines above its own guard.
%
% checkcode passes this file: every line is syntactically perfect.
%
% EXPECTED: use_before_def.py reports cv_seed_mvpa_reg_cov,
%           def line 34 AFTER use line 31.
%
% If a change to the checker makes this file pass, the checker is broken - and
% in a way that reports PASS rather than an error, which is the worse direction.
%
% -------------------------------------------------------------------------
%
% part of the LaBGAScore script-checker positive controls

if ~exist('tuned_seed_mvpa_reg_cov','var') || isempty(tuned_seed_mvpa_reg_cov)
    tuned_seed_mvpa_reg_cov = cv_seed_mvpa_reg_cov;              % <-- USE (line 31)
end

if ~exist('cv_seed_mvpa_reg_cov','var'), cv_seed_mvpa_reg_cov = []; end   % <-- DEF (line 34), too late

if ~isempty(cv_seed_mvpa_reg_cov)
    rng(cv_seed_mvpa_reg_cov, 'twister');
end
