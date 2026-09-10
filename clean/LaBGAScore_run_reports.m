function results = LaBGAScore_run_reports(scriptnames, htmlsavedir, varargin)
% Publish a chain of analysis scripts and FAIL LOUDLY when one of them errors.
%
% *USAGE*
%
% results = LaBGAScore_run_reports(scriptnames, htmlsavedir)
% results = LaBGAScore_run_reports(scriptnames, htmlsavedir, 'artefacts', A)
% results = LaBGAScore_run_reports(..., 'minbytes', 1e6, 'stoponerror', true)
%
% *WHY THIS EXISTS*
%
% MATLAB's publish() CATCHES a script's error, writes it into the html report,
% and then RETURNS NORMALLY. A run that died halfway is therefore
% indistinguishable, from the calling code, from one that succeeded: no
% exception is raised and the exit status is 0. Chains of prep_*/c2* scripts
% have repeatedly appeared to complete while a script had in fact crashed -
% typically on an option variable that only exists in some configurations, or
% after a long computation but before anything was saved.
%
% This function publishes each script and then READS THE REPORT BACK to decide
% whether it actually worked, optionally checking that the artefacts the script
% was supposed to write exist and are of plausible size. It continues past a
% failure rather than aborting the chain, and names every failure at the end.
%
% *HOW FAILURE IS DETECTED*
%
% publish() wraps caught errors in <pre class="codeoutput error"> in the html.
% That MARKUP is what is matched here - deliberately, rather than searching the
% report text for phrases like "Error in" or "Unrecognized function". Those
% phrases occur in ordinary comments and are rendered into the report as normal
% prose, so a text search reports failures for scripts that ran perfectly. This
% was verified: a script whose header comment merely mentions those phrases
% produces them in a <p> block, while its real error appears in the error <pre>.
%
% *INPUTS*
%
%   scriptnames   char or cellstr. Script names WITHOUT the .m extension, in
%                 the order they should run. They must be on the MATLAB path.
%
%   htmlsavedir   directory the reports are written to. Usually the
%                 htmlsavedir defined by the study's a_set_up_paths script.
%
% *OPTIONAL INPUTS*
%
%   'artefacts'   cell array, one entry per script, listing files that script
%                 must have written. Each entry is a char, a cellstr, or []
%                 for "nothing to check". Relative paths are resolved against
%                 'artefactdir'. A script that publishes a clean report but
%                 wrote no results is reported as a failure.
%
%   'artefactdir' base directory for relative artefact paths.
%                 Default: the parent of htmlsavedir, which is the results
%                 directory in the standard LaBGAS layout (results/html).
%
%   'minbytes'    an artefact smaller than this counts as missing.
%                 Default 1000. Catches truncated or placeholder files.
%
%   'stoponerror' true rethrows on the first failure instead of continuing.
%                 Default false - the point of the chain is to learn about
%                 every failure in one pass.
%
%   'publishfcn'  function handle used to publish. Default
%                 @LaBGAScore_prov_publish, which adds the provenance section.
%                 Pass @(nm,dir) publish(nm,'outputDir',dir) for plain publish.
%
% *OUTPUT*
%
%   results       table, one row per script: name, ok, minutes, reportpath,
%                 and reason (empty when ok).
%
% *EXAMPLE*
%
%   scripts = {'proj_secondlevel_m1_s4_prep_2_load_image_data_and_save'
%              'proj_secondlevel_m1_s5_prep_3_calc_univariate_contrasts'
%              'proj_secondlevel_m1_s6_prep_3a_run_second_level_regression'};
%
%   arts    = {'data_objects.mat'
%              'contrast_data_objects.mat'
%              []};
%
%   results = LaBGAScore_run_reports(scripts, htmlsavedir, ...
%                 'artefacts', arts, 'minbytes', 1e6);
%
%   if ~all(results.ok), error('%d report(s) failed', sum(~results.ok)); end
%
% *SEE ALSO*
%
% LaBGAScore_prov_publish, labgascore_run_headless.sh, LaBGAScore_check_display
%
% -------------------------------------------------------------------------
% Lukas Van Oudenhove, KU Leuven, September 2026
% -------------------------------------------------------------------------

% ---------------------------- parse inputs -------------------------------

if ischar(scriptnames) || isstring(scriptnames)
    scriptnames = cellstr(scriptnames);
end
scriptnames = scriptnames(:);
n = numel(scriptnames);

p = inputParser;
p.addParameter('artefacts',   repmat({[]}, n, 1), @iscell);
p.addParameter('artefactdir', '', @(x) ischar(x) || isstring(x));
p.addParameter('minbytes',    1000, @isnumeric);
p.addParameter('stoponerror', false, @islogical);
p.addParameter('publishfcn',  @LaBGAScore_prov_publish, @(x) isa(x,'function_handle'));
p.parse(varargin{:});

artefacts = p.Results.artefacts(:);
if numel(artefacts) ~= n
    error('LaBGAScore_run_reports:artefacts', ...
        '''artefacts'' must have one entry per script (%d given, %d scripts)', ...
        numel(artefacts), n);
end

artefactdir = char(p.Results.artefactdir);
if isempty(artefactdir)
    artefactdir = fileparts(htmlsavedir);   % results/html -> results
end

if ~exist(htmlsavedir, 'dir')
    error('LaBGAScore_run_reports:nohtmldir', ...
        'htmlsavedir does not exist: %s', htmlsavedir);
end

% Fail before doing an hour of work, not after.
missing = scriptnames(cellfun(@(s) isempty(which(s)), scriptnames));
if ~isempty(missing)
    error('LaBGAScore_run_reports:notonpath', ...
        'not on the MATLAB path: %s', strjoin(missing', ', '));
end

% ------------------------------- run -------------------------------------

name       = scriptnames;
ok         = false(n,1);
minutes    = zeros(n,1);
reportpath = repmat({''}, n, 1);
reason     = repmat({''}, n, 1);

fprintf('\n===== publishing %d script(s) to %s =====\n', n, htmlsavedir);

for i = 1:n

    fprintf('\n---------- [%d/%d] %s ----------\n', i, n, name{i});
    t0 = tic;

    try
        h = p.Results.publishfcn(name{i}, htmlsavedir);
        reportpath{i} = h;

        why = local_report_error(h, name{i});

        if isempty(why)
            why = local_artefacts_missing(artefacts{i}, artefactdir, p.Results.minbytes);
        end

        minutes(i) = toc(t0)/60;

        if isempty(why)
            ok(i) = true;
            fprintf('>>> %s OK in %.1f min\n', name{i}, minutes(i));
        else
            reason{i} = why;
            fprintf(2, '>>> %s FAILED after %.1f min: %s\n', name{i}, minutes(i), why);
            if p.Results.stoponerror
                error('LaBGAScore_run_reports:failed', '%s: %s', name{i}, why);
            end
        end

    catch ME
        minutes(i) = toc(t0)/60;
        reason{i}  = ME.message;
        fprintf(2, '>>> %s FAILED after %.1f min: %s\n', name{i}, minutes(i), ME.message);
        if p.Results.stoponerror, rethrow(ME); end
    end

end

results = table(name, ok, minutes, reportpath, reason);

nfail = sum(~results.ok);
fprintf('\n===== %d/%d OK, %d failed, %.1f min total =====\n', ...
    sum(results.ok), n, nfail, sum(results.minutes));
if nfail > 0
    fprintf(2, 'FAILED: %s\n', strjoin(results.name(~results.ok)', ', '));
end

end % main function


% =========================================================================

function why = local_report_error(htmlfile, scriptname)
% Return a description of the error publish() caught into the report, or ''.

why = '';

if isempty(htmlfile) || ~exist(htmlfile, 'file')
    why = 'publish() returned no report file';
    return
end

txt = fileread(htmlfile);

% publish() marks caught errors with an "error" class on the output <pre>.
% Match the markup, NOT the report text: phrases such as "Error in" appear in
% ordinary comments, which publish renders as prose, so a text search flags
% scripts that ran perfectly.
% NB: no \b here. MATLAB's regexp does not treat \b as a word boundary (it
% uses \< and \>), so a pattern like \berror\b silently never matches and
% every failing report is reported as OK. publish() writes exactly
% class="codeoutput error", so a plain substring match is both correct and
% safer than trying to be clever.
blocks = regexp(txt, '<pre class="[^"]*error[^"]*">(.*?)</pre>', 'tokens');

if isempty(blocks), return, end

msg = regexprep(blocks{1}{1}, '<[^>]+>', '');    % strip nested spans
msg = regexprep(msg, '\s+', ' ');
msg = strtrim(local_unescape(msg));

why = sprintf('error in published report (%s): %s', scriptname, ...
    msg(1:min(200, numel(msg))));

end


function why = local_artefacts_missing(want, basedir, minbytes)
% Return a description of missing/too-small artefacts, or ''.

why = '';
if isempty(want), return, end
if ischar(want) || isstring(want), want = cellstr(want); end

bad = {};
for k = 1:numel(want)
    f = char(want{k});
    if ~local_isabs(f), f = fullfile(basedir, f); end
    d = dir(f);
    if isempty(d)
        bad{end+1} = sprintf('%s (missing)', want{k}); %#ok<AGROW>
    elseif d(1).bytes < minbytes
        bad{end+1} = sprintf('%s (%d bytes < %d)', want{k}, d(1).bytes, minbytes); %#ok<AGROW>
    end
end

if ~isempty(bad)
    why = ['report published but expected output not written: ' strjoin(bad, '; ')];
end

end


function tf = local_isabs(f)
tf = ~isempty(f) && (f(1) == filesep || (ispc && numel(f) > 1 && f(2) == ':'));
end


function s = local_unescape(s)
s = strrep(s, '&lt;',  '<');
s = strrep(s, '&gt;',  '>');
s = strrep(s, '&quot;','"');
s = strrep(s, '&#39;', '''');
s = strrep(s, '&amp;', '&');
end
