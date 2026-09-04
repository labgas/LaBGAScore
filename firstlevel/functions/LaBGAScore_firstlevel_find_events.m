function [evf, how] = LaBGAScore_firstlevel_find_events(BIDSdir, sub, stem)
% LaBGAScore_firstlevel_find_events  Resolve one run's events.tsv by BIDS inheritance.
%
%   evf = LaBGAScore_firstlevel_find_events(BIDSdir, sub, stem)
%   [evf, how] = LaBGAScore_firstlevel_find_events(BIDSdir, sub, stem)
%
%   Returns the full path to the events.tsv that applies to one functional run,
%   searching from the most specific level of the BIDS hierarchy upwards, or ''
%   when nothing applies.
%
%   Three routes are tried, most trustworthy first:
%
%     exact        <BIDSdir>/<sub>/func/<stem>_events.tsv, i.e. a run-specific
%                  file sitting next to the functional image
%     inheritance  the BIDS inheritance principle: a file higher up the tree
%                  applies to everything below it, provided every entity in its
%                  own name matches this run. A block design with fixed timing
%                  is often stored once as <BIDSdir>/task-<label>_events.tsv
%                  rather than copied per subject and run. Session, subject and
%                  dataset levels are searched in that order, and the most
%                  specific applicable candidate wins
%     run entity   a single events file in the subject's func dir carrying the
%                  right run, used when the task label is spelled differently in
%                  BIDS and in derivatives (e.g. task-emosex_movies against
%                  task-emosex), which no entity match would survive
%
%   INPUTS
%     BIDSdir  root of the BIDS dataset
%     sub      subject directory name, e.g. 'sub-01'
%     stem     filename stem of the run, without suffix - e.g.
%              'sub-01_task-MIST_run-1', as built by the calling script from the
%              fMRIPrep confounds filename
%
%   OUTPUTS
%     evf  full path to the events file, or '' when none applies
%     how  which route matched: 'exact', 'inheritance', 'run entity only', or ''
%
%   Files whose name contains 'noninterest' are ignored throughout, since
%   LaBGAS stores nuisance events alongside the events of interest.
%
%   This is the resolution used by LaBGAScore_firstlevel_task_motion_diagnostics,
%   factored out so the fitting scripts resolve events the same way.
%
%   See also LaBGAScore_firstlevel_s2_fit_model,
%   LaBGAScore_firstlevel_s2a_fit_model_multisess_multitask.

% AUTHOR
% Lukas Van Oudenhove, KU Leuven
%
% ..

evf = ''; how = '';

% 1. exact: a run-specific file next to the functional image
exact = fullfile(BIDSdir, sub, 'func', [stem '_events.tsv']);
if isfile(exact), evf = exact; how = 'exact'; return, end

want = local_entities(stem);

% 2. inheritance: search from the most specific level upwards
levels = { fullfile(BIDSdir, sub, 'func'), fullfile(BIDSdir, sub), BIDSdir };
if isfield(want, 'ses')
    levels = [ { fullfile(BIDSdir, sub, ['ses-' want.ses], 'func'), ...
                 fullfile(BIDSdir, sub, ['ses-' want.ses]) }, levels ];
end

for L = 1:numel(levels)
    if ~isfolder(levels{L}), continue, end
    cand = dir(fullfile(levels{L}, '*_events.tsv'));
    if isempty(cand), continue, end
    cand = cand(~[cand.isdir] & ~contains({cand.name}, 'noninterest'));
    best = ''; bestn = -1;
    for c = 1:numel(cand)
        have = local_entities(cand(c).name);
        keys = fieldnames(have);
        ok = true;
        for k = 1:numel(keys)
            if ~isfield(want, keys{k}) || ~strcmp(want.(keys{k}), have.(keys{k}))
                ok = false; break
            end
        end
        % a candidate with more matching entities is the more specific one
        if ok && numel(keys) > bestn
            best = fullfile(cand(c).folder, cand(c).name); bestn = numel(keys);
        end
    end
    if ~isempty(best), evf = best; how = 'inheritance'; return, end
end

% 3. last resort: right run, whatever the task label
runtok = regexp(stem, 'run-[0-9]+', 'match', 'once');
cand = dir(fullfile(BIDSdir, sub, 'func', '*_events.tsv'));
if isempty(cand), return, end
cand = cand(~contains({cand.name}, 'noninterest'));
if ~isempty(runtok), cand = cand(contains({cand.name}, runtok)); end
if isscalar(cand)
    evf = fullfile(cand(1).folder, cand(1).name); how = 'run entity only';
end

end % function


function e = local_entities(name)
% BIDS key-value entities from a filename, as a struct
    e = struct();
    parts = strsplit(erase(name, '.tsv'), '_');
    for i = 1:numel(parts)
        kv = strsplit(parts{i}, '-');
        if numel(kv) == 2, e.(matlab.lang.makeValidName(kv{1})) = kv{2}; end
    end
end
