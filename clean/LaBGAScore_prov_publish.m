function htmlfile = LaBGAScore_prov_publish(scriptname, htmlsavedir, varargin)
%
%
% *USAGE*
%
% Drop-in replacement for the publish() call documented in the header of
% every LaBGAS second-level script. Where you would type
%
%   publish('cfs_secondlevel_m14_s7_c2a_second_level_regression', ...
%           'outputDir', htmlsavedir)
%
% type instead
%
%   LaBGAScore_prov_publish('cfs_secondlevel_m14_s7_c2a_second_level_regression', ...
%                           htmlsavedir)
%
% and the html report gains a Provenance section recording the commit of
% every dependency the script reaches, which of them carried uncommitted
% local changes, and the MATLAB/SPM versions in play. A machine-readable
% copy is written to the results notes directory at the same time.
%
% Deliberately a wrapper rather than an edit to the ~40 script templates:
% the templates live partly in this repo and partly in the LaBGAS fork of
% CANlab_help_examples, and each study has its own renamed copies that
% would never pick up such an edit. Changing what you type changes nothing
% about the scripts themselves.
%
% publish() evaluates code in the base workspace, so wrapping it in a
% function does not disturb the workflow's reliance on variables inherited
% from a_set_up_paths_always_run_first and the prep scripts.
%
%
% *OPTIONS*
%
% * scriptname      script to publish, name or full path
% * htmlsavedir     output directory. Default: htmlsavedir from the base
%                   workspace.
% * 'savedir'       where the machine-readable record goes. Default:
%                   notesdir from the base workspace.
% * 'showCode'      default true, as in the existing LaBGAS publish calls
% * 'maxWidth'      NOT set by default. publish() honours this by resizing
%                   the .png on disk, so detail is permanently lost rather
%                   than merely displayed smaller: a 1600x900 figure
%                   published with maxWidth 1600 / maxHeight 800 lands as
%                   1422x800. Leaving it unset reproduces exactly what the
%                   bare publish() call documented in each script header
%                   already does. Set it only when you deliberately want
%                   smaller files.
% * 'maxHeight'     as above, not set by default
% * (the dimensions of every figure publish() writes are recorded in
%    PROV.figures and appended to the .mat record, alongside the screen
%    geometry, so the same script run on another machine can be compared)
% * 'responsive'    default true. Injects a viewport tag and a small
%                   stylesheet so figures scale down on narrow screens and
%                   wide console output scrolls in its own box instead of
%                   stretching the page. Set false for byte-identical
%                   publish() output.
% * 'snapshot'      set false to publish without recording provenance
% * any other name-value pair is passed through to publish()
%
%
% *OUTPUT*
%
% * htmlfile        full path to the html report produced
%
%
% *NOTES*
%
% The provenance table is injected into the html AFTER publish() returns.
% That is safe here because the file has just been created and is not yet
% under datalad control; the retrospective tool deliberately does NOT do
% this to existing reports, and writes sidecar files instead.
%
%
% *DISPLAY SETTINGS*
%
% publish() captures figures from the screen, so your X2go window size and
% display DPI decide how figures in the report come out. The header of every
% report produced here records both, plus the resulting figure dimensions,
% and flags figures whose size was set by the display rather than by the
% script. Run LaBGAScore_check_display to check your session, and see
% "Set up your X2go display for publishing figures" in
% LaBGAS_fMRI_analysis_workflow.md for the recommended settings per screen.
%
%
% *SEE ALSO*
%
% LaBGAScore_prov_snapshot, LaBGAScore_prov_resolve_retrospective,
% LaBGAScore_check_display, plugin_set_figure_size (CANlab_help_examples)
%
% -------------------------------------------------------------------------
%
% modified by: Lukas Van Oudenhove & Claude Opus 5
%
% date:   KU Leuven, August, 2026
%
% -------------------------------------------------------------------------
%
% LaBGAScore_prov_publish.m         v1.0
%
% last modified: 2026/08/31
%
%
%% PARSE OPTIONS
% -------------------------------------------------------------------------

p = inputParser;
p.KeepUnmatched = true;
p.addParameter('savedir', '__default__', @(x) ischar(x) || isstring(x));
p.addParameter('showCode', true, @islogical);
p.addParameter('maxWidth', [], @isnumeric);
p.addParameter('maxHeight', [], @isnumeric);
p.addParameter('responsive', true, @islogical);
p.addParameter('snapshot', true, @islogical);
p.parse(varargin{:});
opt = p.Results;

if nargin < 2 || isempty(htmlsavedir)
    try
        htmlsavedir = evalin('base','htmlsavedir');
    catch
        error(['htmlsavedir not given and not found in the base workspace. ' ...
               'Run a_set_up_paths_always_run_first first, or pass it explicitly.']);
    end
end

htmlsavedir = char(htmlsavedir);
if ~isfolder(htmlsavedir), mkdir(htmlsavedir); end

[~, scriptbase] = fileparts(char(scriptname));


%% TAKE THE SNAPSHOT BEFORE RUNNING
% -------------------------------------------------------------------------
% before, not after: the script may itself modify files, and what we want
% on record is the state the script ran against

PROV = [];

if opt.snapshot
    try
        PROV = LaBGAScore_prov_snapshot('scriptname', char(scriptname), ...
            'savedir', opt.savedir, 'print', false);
    catch ME
        warning('LaBGAScore_prov_publish:snapshot', ...
            'provenance snapshot failed (%s); publishing without it', ME.message);
    end
end


%% PUBLISH
% -------------------------------------------------------------------------

pubopts = struct('useNewFigure', false, ...
                 'format', 'html', ...
                 'outputDir', htmlsavedir, ...
                 'showCode', opt.showCode);

% Only pass the size caps when they are asked for. They are NOT set by
% default, because publish() applies them by resizing the .png on disk - the
% detail is gone, not merely displayed smaller. A 1600x900 figure published
% with maxWidth 1600 / maxHeight 800 lands as 1422x800. Brain montages here
% are 1707x913 natively, which is where that detail matters most. Fitting
% small screens is a display problem and is solved by the responsive CSS
% below, not by throwing away pixels for every reader.
if ~isempty(opt.maxWidth),  pubopts.maxWidth  = opt.maxWidth;  end
if ~isempty(opt.maxHeight), pubopts.maxHeight = opt.maxHeight; end

extra = fieldnames(p.Unmatched);
for k = 1:numel(extra)
    pubopts.(extra{k}) = p.Unmatched.(extra{k});
end

htmlfile = publish(char(scriptname), pubopts);


%% RECORD WHAT THE FIGURES ACTUALLY CAME OUT AS
% -------------------------------------------------------------------------
% Figure capture depends on the display: publish() grabs what is on screen,
% so the same script run from a different machine produces different pixel
% dimensions. Recording the result next to the screen geometry turns that
% from a puzzling difference between two reports into a measured one.

if ~isempty(PROV)
    try
        PROV.figures = local_figure_summary(htmlsavedir, scriptbase, PROV);
        if ~isempty(PROV.files_written)
            matfile = char(PROV.files_written(1));
            if isfile(matfile)
                save(matfile, 'PROV', '-append');
            end
        end
    catch ME
        warning('LaBGAScore_prov_publish:figures', ...
            'could not record figure dimensions (%s)', ME.message);
    end
end


%% INJECT THE PROVENANCE SECTION INTO THE REPORT
% -------------------------------------------------------------------------

if isfile(htmlfile) && opt.responsive
    try
        local_make_responsive(htmlfile);
    catch ME
        warning('LaBGAScore_prov_publish:responsive', ...
            'could not add the responsive stylesheet to %s (%s)', htmlfile, ME.message);
    end
end

if ~isempty(PROV) && isfile(htmlfile)
    try
        local_inject(htmlfile, PROV, scriptbase);
    catch ME
        warning('LaBGAScore_prov_publish:inject', ...
            'could not add the provenance section to %s (%s)', htmlfile, ME.message);
    end
end

fprintf('\nPublished %s\n  report     %s\n', scriptbase, htmlfile);
if ~isempty(PROV) && isfield(PROV,'figures') && ~isempty(PROV.figures)
    fprintf('  figures    %d, up to %dx%d px\n', height(PROV.figures), ...
        max(PROV.figures.Width), max(PROV.figures.Height));
    if any(PROV.figures.ScreenLimited)
        fprintf(['  NOTE       %d figure(s) are as wide as this display, so their size was\n' ...
                 '             set by the screen rather than by the script. On a different\n' ...
                 '             machine this report will look different. See\n' ...
                 '             plugin_set_figure_size in CANlab_help_examples.\n'], ...
                 sum(PROV.figures.ScreenLimited));
    end
end
if ~isempty(PROV) && ~isempty(PROV.files_written)
    fprintf('  provenance %s\n', PROV.files_written(1));
end
fprintf('\n');

end


%% LOCAL FUNCTIONS
% -------------------------------------------------------------------------

function F = local_figure_summary(htmlsavedir, scriptbase, PROV)
% Dimensions of every figure publish() wrote for this script, and whether
% any of them looks limited by the size of the display rather than by the
% size that was asked for.

F = table(strings(0,1), zeros(0,1), zeros(0,1), false(0,1), ...
    'VariableNames', {'File','Width','Height','ScreenLimited'});

d = dir(fullfile(htmlsavedir, [scriptbase '_*.png']));
d = d(~[d.isdir]);

if isempty(d), return, end

screenw = NaN;
if isfield(PROV,'env') && isfield(PROV.env,'screen_size')
    wh = sscanf(char(PROV.env.screen_size), '%dx%d');
    if numel(wh) == 2, screenw = wh(1); end
end

for k = 1:numel(d)
    fp = fullfile(d(k).folder, d(k).name);
    try
        info = imfinfo(fp);
    catch
        continue
    end
    % within a few percent of the full screen width means the display, not
    % the requested size, decided how big this figure is
    limited = ~isnan(screenw) && info(1).Width >= 0.97 * screenw;
    F = [F; {string(d(k).name), info(1).Width, info(1).Height, limited}]; %#ok<AGROW>
end

end


function local_make_responsive(htmlfile)
% Make a published report usable on any screen.
%
% publish() emits <img> tags with no width attributes, a stylesheet with no
% max-width on images and no width constraint on .content, and no viewport
% meta tag. A 1707px brain montage therefore overflows a 1366px laptop and
% scrolls the whole page sideways. These overrides are injected after
% publish's own <style>, so they win, and they change nothing on a screen
% wide enough to show the figure at native size.

html = fileread(htmlfile);

if contains(html, 'labgascore-responsive')
    return                                  % already done
end

nl = newline;

css = [ ...
    '<meta name="viewport" content="width=device-width, initial-scale=1">' nl ...
    '<style type="text/css" id="labgascore-responsive">' nl ...
    '/* Added by LaBGAScore_prov_publish. Figures keep their full resolution' nl ...
    '   on disk; this only stops them overflowing narrow screens. */' nl ...
    'img { max-width: 100%; height: auto; }' nl ...
    '/* long console output scrolls inside its own box instead of' nl ...
    '   stretching the page */' nl ...
    'pre.codeinput, pre.codeoutput, pre.error { overflow-x: auto; }' nl ...
    '/* wide stats tables likewise */' nl ...
    '.labgas-scroll { overflow-x: auto; }' nl ...
    '</style>' nl];

idx = strfind(html, '</head>');

if isempty(idx)
    return
end

html = [html(1:idx(1)-1) css html(idx(1):end)];

fid = fopen(htmlfile, 'w');
if fid < 0, error('cannot write %s', htmlfile); end
fwrite(fid, html);
fclose(fid);

end


function local_inject(htmlfile, PROV, scriptbase)

html = fileread(htmlfile);

block = local_html_block(PROV, scriptbase);

% put it at the end of the body, ahead of the source listing comment that
% publish appends after </body>
idx = strfind(html, '</body>');

if isempty(idx)
    html = [html block];
else
    html = [html(1:idx(end)-1) block html(idx(end):end)];
end

fid = fopen(htmlfile, 'w');
if fid < 0, error('cannot write %s', htmlfile); end
fwrite(fid, html);
fclose(fid);

end


function s = local_html_block(PROV, scriptbase)

nl = newline;

used = PROV.deps(PROV.deps.UsedByScript,:);
if isempty(used), used = PROV.deps; end

% only repositories this script calls directly, or that carry a relevant
% uncommitted change. A repo reached at depth 6 through a chain of ambiguous
% names is not a dependency anyone would recognise as one.
hasdepmap = ismember('MinDepth', used.Properties.VariableNames) && any(used.NFilesUsed > 0);
if ismember('MinDepth', used.Properties.VariableNames)
    nindirect = sum(~(used.MinDepth <= 2 | used.DirtyRelevant));
    used = used(used.MinDepth <= 2 | used.DirtyRelevant,:);
else
    nindirect = 0;
end

s = [nl '<h2 id="provenance">Provenance</h2>' nl];
s = [s '<p style="font-size:90%">Recorded by LaBGAScore_prov_publish at run time. ' ...
     'Commit hashes identify the state of each dependency when this report was produced.</p>' nl];

s = [s '<p style="font-size:90%"><b>' char(scriptbase) '</b> &middot; ' ...
     char(PROV.env.timestamp) ' &middot; ' char(PROV.env.hostname) ...
     ' &middot; MATLAB ' char(PROV.env.matlab_release)];
if isfield(PROV.env,'screen_size')
    % figure sizes in this report depend on this screen - see
    % plugin_set_figure_size in CANlab_help_examples
    s = [s ' &middot; screen ' char(PROV.env.screen_size) ...
         ' @ ' num2str(PROV.env.screen_dpi) ' DPI'];
end
if isfield(PROV,'figures') && ~isempty(PROV.figures)
    s = [s ' &middot; figures up to ' ...
         sprintf('%dx%d', max(PROV.figures.Width), max(PROV.figures.Height))];
    if any(PROV.figures.ScreenLimited)
        s = [s ' <span style="color:#b00">(display-limited)</span>'];
    end
end
for k = 1:height(PROV.nonrepo)
    s = [s ' &middot; ' char(PROV.nonrepo.Name(k)) ' ' char(PROV.nonrepo.Version(k))]; %#ok<AGROW>
end
s = [s '</p>' nl];

s = [s '<div class="labgas-scroll">' nl];
s = [s '<table cellpadding="4" style="border-collapse:collapse;font-size:90%">' nl];
s = [s '<tr style="border-bottom:1px solid #999"><th align="left">Repository</th>' ...
     '<th align="left">Commit</th><th align="left">Branch</th>' ...
     '<th align="right">Files used</th><th align="left">State</th></tr>' nl];

for k = 1:height(used)
    if used.DirtyRelevant(k)
        state = ['<span style="color:#b00"><b>uncommitted changes in files this script uses:</b> ' ...
                 char(strrep(used.RelevantModifiedFiles(k),'; ','<br>')) '</span>'];
    elseif used.Dirty(k) && ~hasdepmap
        state = sprintf('%d modified file(s); relevance to this script unknown', used.NModified(k));
    elseif used.Dirty(k)
        state = sprintf('clean for this script (%d modified elsewhere in repo)', used.NModified(k));
    else
        state = 'clean';
    end
    s = [s sprintf(['<tr style="border-bottom:1px solid #ddd"><td>%s</td><td><code>%s</code></td>' ...
                    '<td>%s</td><td align="right">%d</td><td>%s</td></tr>' nl], ...
        char(used.Repo(k)), char(used.CommitShort(k)), char(used.Branch(k)), ...
        used.NFilesUsed(k), state)]; %#ok<AGROW>
end

s = [s '</table>' nl '</div>' nl];

if nindirect > 0
    s = [s sprintf(['<p style="font-size:85%%;color:#666">%d further repositories are ' ...
         'reachable only indirectly (depth &gt;2); the full list is in the .mat ' ...
         'provenance record.</p>%s'], nindirect, nl)];
end

if any(used.DirtyRelevant)
    s = [s '<p style="font-size:90%;color:#b00">The commit hashes above do not fully ' ...
         'describe this run: the flagged files carried uncommitted local edits. ' ...
         'Their exact content at run time is stored in the .mat provenance record ' ...
         'alongside this report.</p>' nl];
end

end
