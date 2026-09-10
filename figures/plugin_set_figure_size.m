function [fh, actual_size] = plugin_set_figure_size(varargin)
% Size a figure for capture into a publish() HTML report.
%
% Lives in LaBGAScore/figures/ and is used by BOTH the first-level scripts in
% LaBGAScore and the second-level templates in the LaBGAS fork of
% CANlab_help_examples. It started life in that fork, as one of its plugin_*
% helpers, and moved here when first level adopted it: a general-purpose
% figure utility belongs in the core repo, and the fork already depends on
% LaBGAScore for far more than this.
%
% The plugin_ prefix is kept deliberately, against LaBGAScore's usual
% LaBGAScore_* convention. MATLAB resolves it by name, so leaving the name
% alone moved the file without touching any of its 144 call sites; renaming it
% would have meant editing all of them for nothing.
%
% Sizes the current figure (gcf) for capture into a publish() HTML report,
% choosing the largest canvas that the current session's display can
% actually capture. Replaces set(gcf,'WindowState','maximized').
%
% WHY NOT 'maximized': maximizing ties figure pixel dimensions to whichever
% client screen (e.g. an X2go session) happened to be connected when the
% script ran. All figures in these scripts use MATLAB's default font sizes,
% which are specified in POINTS (1 pt = 1/72 inch, a fixed physical unit).
% A canvas whose physical size changes per session therefore renders text
% inconsistently too large or too small depending on who ran the script.
%
% WHY INCHES, NOT PIXELS: anchoring the canvas in inches - a physical unit
% related to points by an exact, DPI-independent conversion - keeps the
% font-to-canvas ratio constant across sessions whose ScreenPixelsPerInch
% differs. A pixel-fixed canvas would not: point-based text would render at
% a different pixel size on each session while the canvas stayed put. Do
% not change the primary specification back to pixels.
%
% WHY IT NOW FITS TO THE SCREEN: publish() captures what is on screen. A
% figure larger than the display - or positioned partly off it - is
% captured at display size instead, AND at a different aspect ratio than
% requested, silently. get(fh,'Position') still reports the size you asked
% for, so nothing warns you. Measured on the LaBGAS server (1718x1360 at
% 133 DPI): the previous fixed 16x10 inch default needs 2128x1330 px, does
% not fit, and captured 1718x1254 - aspect 1.37 instead of 1.60, i.e. the
% same result as 'maximized', which is what this function exists to avoid.
% The failure appears on high-DPI sessions, which is precisely the case it
% was written for.
%
% WHY THE DEFAULT IS 12 x 7.5 INCHES, NOT 16 x 10: 16 x 10 in is not
% reachable on any lab laptop. It needs 16*DPI x 10*DPI pixels of window,
% so a 1366x768 client would have to run at 72 DPI and a 1600x900 client at
% 84 DPI - both impractically small to work in. 12 x 7.5 in (same 16:10
% aspect) is reachable on every screen in the lab at 96 DPI, and leaves
% headroom up to ~135 DPI on a 1920x1080 client and ~140 DPI on half of a
% 3440x1440 ultrawide. A default nobody can actually achieve guarantees the
% inconsistency this function exists to prevent. Run
% LaBGAScore_check_display (LaBGAScore/clean) to see what your own session
% can do, and pass 'width'/'height' explicitly if you want something else.
%
% So the requested size is now treated as an upper bound. If it does not
% fit the display, both dimensions are scaled down by the same factor, so
% the ASPECT RATIO IS ALWAYS HONOURED and only the absolute size gives way.
% The font-to-canvas ratio is then constant across every session whose
% display can hold the requested size, and degrades gracefully (text
% relatively larger, which keeps a smaller capture legible) below that.
% The figure is also repositioned fully on-screen, since a window hanging
% off the edge is clamped at capture no matter how it was sized.
%
% Call this immediately before drawnow/snapnow (and before
% plugin_save_figure, if used), after all plotting/display calls for the
% figure are complete.
%
% USAGE:
% plugin_set_figure_size()                         % default, fitted to screen
% plugin_set_figure_size('nrows', nrows)           % multi-row montage (canlab_results_fmridisplay 'multirow')
% plugin_set_figure_size('width', w, 'height', h)  % explicit upper bound
% plugin_set_figure_size('fig', fh)                % a figure a drawing call opened
% plugin_set_figure_size('fig', fh, 'minpanel', [1.2 1.0])   % grid of many panels
% plugin_set_figure_size('titlescale', 0.5)        % dense layout, smaller titles
%
% OPTIONAL NAME-VALUE ARGUMENTS:
% 'width'   maximum figure width in inches (default 12)
% 'height'  maximum figure height in inches (default 7.5 if 'nrows' not
%           given; keeps the 16:10 aspect of the previous 16x10 default)
% 'nrows'   number of montage rows passed to canlab_results_fmridisplay's
%           'multirow' option (e.g. num_effects). Does NOT scale height:
%           canlab_results_fmridisplay's own 'multirow' code allocates each
%           row a FIXED normalized-coordinate band and sizes the whole
%           figure the same way regardless of how many rows it holds (1-4
%           per figure, extra rows spill into new figures) - so a canvas
%           shrunk for fewer rows starves every row of physical space
%           rather than saving any. Confirmed via real-project testing: an
%           earlier version that shrank height for low nrows produced
%           montages with the title clipped and slices squeezed into a
%           sliver at the top. 'nrows' is accepted (so multirow call sites
%           can still document their row count) but currently only uses the
%           same default height as any other figure.
% 'margin'  [horizontal vertical] fraction of the screen to leave free for
%           window decorations and panels (default [0.02 0.06])
% 'minsize' warn if fitting forces the canvas below this width in inches
%           (default 7). A capture much smaller than this makes montage
%           text hard to read in the report.
% 'verbose' print a line when the requested size had to be reduced
%           (default true). Reported only ONCE per session per distinct
%           request/display combination, since a published report calls
%           this once per figure. Silence it entirely with false.
% 'fig'     handle(s) of the figure(s) to size (default gcf). Use this for
%           drawing calls that OPEN THEIR OWN figure - canlab_results_fmridisplay
%           with 'multirow', @region/montage, and histogram(...,'byimage') all
%           do - where gcf is no longer the figure you want to size. Pass a
%           vector to size several at once.
% 'titlescale'
%           factor applied to the font size of every title in the figure
%           (default 2/3). Titles scale with the canvas, so a canvas grown to
%           fit many panels would otherwise carry enormous titles. Pass a
%           smaller value for dense layouts: the carpet-plot panels use 0.5,
%           whose 15.4 pt titles overlap even at the 2/3 default. Each title is
%           marked once it has been scaled, so titles added after this call are
%           still scaled on a later call and none is scaled twice.
% 'keepaspect'
%           true preserves the figure's CURRENT aspect ratio instead of the
%           requested one, scaling it to fit (default false). For figures whose
%           layout is meaningful and unusual - a wide 1x3 strip, say - which
%           would be distorted by being forced to 16:10.
% 'minpanel'
%           [w h] minimum size in inches for a single panel of a grid figure.
%           Grows the canvas until the median panel is at least this big,
%           measuring the ACTUAL layout rather than assuming one, so a
%           per-subject density grid stays legible with 64 or 158 subjects
%           instead of being squeezed into the default canvas. Bounded: the
%           canvas is capped at 20 x 30 inches and a request whose aspect ratio
%           would exceed 3:1 is refused with an explanation, since a layout with
%           that many panels in one row cannot be made readable by resizing and
%           should be split across figures instead.
%
% HEADLESS: when there is no display (feature('ShowFigureWindows') == 0),
% publish() PRINTS figures rather than capturing them from the screen, so the
% fit-to-screen logic above does not apply and figure size is not limited by the
% 1024x768 / 72 dpi virtual screen headless MATLAB reports. The screen bound is
% therefore lifted headless, but ONLY when 'minpanel' is given - deliberately, so
% that every other figure comes out at exactly the same size as before.
%
% OUTPUT:
% fh            handle of the resized figure
% actual_size   [width height] in inches actually applied
%
% Example:
% figure; plot(1:10);
% plugin_set_figure_size();
% drawnow, snapnow;
%
% See also: LaBGAScore_prov_publish (LaBGAScore/clean), which records the
% session's screen size and DPI in the published report, so figure
% differences between machines are diagnosable after the fact.

p = inputParser;
addParameter(p, 'width', 12, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'height', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
addParameter(p, 'nrows', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
addParameter(p, 'margin', [0.02 0.06], @(x) isnumeric(x) && numel(x) == 2 && all(x >= 0 & x < 0.5));
addParameter(p, 'minsize', 7, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'verbose', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'fig', [], @(x) isempty(x) || all(isgraphics(x)));
addParameter(p, 'keepaspect', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'titlescale', 2/3, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'minpanel', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==2 && all(x>0)));
parse(p, varargin{:});

width = p.Results.width;
height = p.Results.height;
margin = p.Results.margin;
titlescale = p.Results.titlescale;

if isempty(height)
    height = 7.5;                      % keeps the 16:10 aspect of the old 16x10 default
end

fh = p.Results.fig;
if isempty(fh)
    fh = gcf;                          % default: the current figure
end

% Several CANlab drawing calls open MORE THAN ONE figure - plot(fmri_data), for
% instance, opens both 'canlab_orthviews' and 'fmri data matrix'. Sizing only gcf
% leaves the others at whatever size they were created with. Accept a vector and
% handle each in turn.
if numel(fh) > 1
    args = varargin;
    for i_fh = 1:numel(fh)
        a = args;
        k = find(strcmpi(a, 'fig'));
        if isempty(k), a = [a, {'fig', fh(i_fh)}]; else, a{k+1} = fh(i_fh); end %#ok<AGROW>
        plugin_set_figure_size(a{:});
    end
    return
end

% WindowState 'maximized' silently overrides any Position set while it is
% active, so it must be cleared first.
set(fh, 'WindowState', 'normal');


%% FIT THE REQUEST TO WHAT THIS DISPLAY CAN CAPTURE
% -------------------------------------------------------------------------

screen_px = get(0, 'ScreenSize');           % [1 1 width height], pixels
dpi = get(0, 'ScreenPixelsPerInch');

% HEADLESS: no clamping. The fitting below exists because publish() run from a
% DESKTOP session captures what is on screen, so a figure bigger than the display
% is silently captured at display size and at the wrong aspect. Under -nodisplay
% there is no screen capture: publish() prints the figure, and the PNG comes out
% at exactly the requested inches x 72 dpi. Measured: 12x7.5 in -> 864x540,
% 20x15 -> 1440x1080, 26x20 -> 1872x1440, all well beyond the 1024x768 virtual
% screen. Clamping headless therefore throws away resolution for no reason - and
% it is precisely the multi-panel figures (one density plot per subject) that
% need to grow past it.
% Deliberately narrow: only figures that ASK to grow (minpanel) are un-clamped.
% Every other figure keeps the exact screen-fitted size it had before, so this
% change cannot move anything that is already correct in the published reports.
is_headless = ~feature('ShowFigureWindows');

if is_headless && ~isempty(p.Results.minpanel)
    % A generous but FINITE bound. Infinity would make 'keepaspect' - which grows
    % a figure by min(usable/current) - grow it without limit.
    usable_w_in = 40;
    usable_h_in = 40;
else
    usable_w_in = screen_px(3) * (1 - margin(1)) / dpi;
    usable_h_in = screen_px(4) * (1 - margin(2)) / dpi;
end

% 'keepaspect': enlarge the figure at ITS OWN aspect ratio rather than forcing the
% default 16:10. Wide, short figures such as canlab_orthviews (819 x 292 px, aspect
% 2.8) are distorted by a 12 x 7.5 in canvas; they want the same shape, bigger.
if p.Results.keepaspect
    set(fh, 'Units', 'inches');
    cur = get(fh, 'Position');
    if cur(3) > 0 && cur(4) > 0
        width  = cur(3);
        height = cur(4);
        grow = min(usable_w_in / width, usable_h_in / height);

        if grow > 1
            width  = width  * grow;
            height = height * grow;
        end
    end
end

% 'minpanel': grow the canvas so that every panel of a multi-panel figure gets at
% least [w h] inches. Read from the axes actually present rather than from any
% assumption about the layout, so it adapts to however many subjects the caller
% happened to plot - a per-subject density plot grid with 158 subjects needs a
% far taller canvas than one with 20, and nobody should have to hand-tune that.
if ~isempty(p.Results.minpanel)

    ax_mp = findobj(fh, 'Type', 'axes');

    if ~isempty(ax_mp)

        pos_mp = get(ax_mp, 'Position');
        if iscell(pos_mp), pos_mp = cell2mat(pos_mp); end

        med_w = median(pos_mp(:,3));    % panel width  as a fraction of the figure
        med_h = median(pos_mp(:,4));    % panel height as a fraction of the figure

        want_w = width;  if med_w > 0, want_w = max(width,  p.Results.minpanel(1) / med_w); end
        want_h = height; if med_h > 0, want_h = max(height, p.Results.minpanel(2) / med_h); end

        % SANITY BOUNDS. A layout that is one long row (many columns, one row)
        % asks for a canvas hundreds of inches wide, which just pins the figure
        % to the headless bound and produces an unreadable ribbon - a real
        % 2880 x 205 px figure came out of prep_3 this way. Cap each dimension,
        % and refuse a wildly non-rectangular canvas outright rather than emit
        % something unusable.
        max_w = 20; max_h = 30; max_aspect = 3;

        want_w = min(want_w, max_w);
        want_h = min(want_h, max_h);

        if want_w / want_h > max_aspect || want_h / want_w > max_aspect
            if p.Results.verbose
                fprintf(['plugin_set_figure_size: minpanel would need a %.0f x %.0f in canvas ' ...
                         '(aspect %.1f) for this layout, which is not usable; leaving the ' ...
                         'default size. The panels are too many for one figure - consider ' ...
                         'plotting fewer images per figure.\n'], want_w, want_h, max(want_w/want_h, want_h/want_w));
            end
        else
            width  = want_w;
            height = want_h;
        end

    end

end

% one scale factor for both dimensions, so the aspect ratio survives
scale = min([1, usable_w_in / width, usable_h_in / height]);

actual_size = [width height] * scale;

% Report at most once per session per distinct situation. A published
% report calls this once per figure - often eight or more times - and the
% same notice repeated down the page is noise, not information.
persistent announced
if isempty(announced), announced = {}; end

situation = sprintf('%g_%g_%d_%d_%g', width, height, screen_px(3), screen_px(4), dpi);
firsttime = ~ismember(situation, announced);

if firsttime
    announced{end+1} = situation;
end

if p.Results.verbose && scale < 1 && firsttime
    fprintf(['plugin_set_figure_size: %.3g x %.3g in does not fit this display ' ...
             '(%d x %d px at %g DPI); using %.3g x %.3g in instead, aspect ratio ' ...
             'preserved. Reported once per session.\n'], width, height, ...
             screen_px(3), screen_px(4), dpi, actual_size(1), actual_size(2));
end

if p.Results.verbose && actual_size(1) < p.Results.minsize && firsttime
    warning('plugin_set_figure_size:smallCanvas', ...
        ['this display only allows a %.3g in wide figure, below the %.3g in ' ...
         'guideline; montage text may be hard to read in the published report. ' ...
         'Consider running from a session with a larger or lower-DPI display.'], ...
        actual_size(1), p.Results.minsize);
end


%% APPLY, KEEPING THE WHOLE WINDOW ON SCREEN
% -------------------------------------------------------------------------
% A window that extends past the screen edge is clamped at capture just as a
% too-large one is, so position matters as much as size. Anchor near the
% bottom-left, which is safe for any size that fits.

set(fh, 'Units', 'inches');

left = margin(1) * screen_px(3) / dpi / 2;
bottom = margin(2) * screen_px(4) / dpi / 2;

set(fh, 'Position', [left, bottom, actual_size(1), actual_size(2)]);


%% SCALE MONTAGE TITLE FONTS
% -------------------------------------------------------------------------
% CanlabCore's title_montage hardcodes FontSize 18
% (@fmridisplay/title_montage), a size tuned for a maximized window. On the
% 12 x 7.5 in canvas this function produces, 18 pt titles are too heavy - and
% on 'regioncenters' montages they are worse still, because @region/montage
% calls title_montage once per region, putting an 18 pt title over each small
% per-region axis. Scale them here rather than in CanlabCore, which is shared
% and deliberately left untouched.
%
% Scaled once per figure: a second call on the same figure would compound the
% reduction, and some scripts size a figure more than once.

if titlescale ~= 1

    ax = findobj(fh, 'Type', 'axes');

    for i = 1:numel(ax)

        th = get(ax(i), 'Title');

        if ~isempty(th) && all(isgraphics(th)) && ~isempty(get(th, 'String'))

            % Mark each TITLE, not the figure. A figure-level flag stops titles
            % that are added AFTER the first sizing call from ever being scaled -
            % which is what happened when one montage block drew into another's
            % figure: the figure was already flagged, so the new 18 pt title was
            % skipped. Per-title marking still prevents double-shrinking.
            if ~isappdata(th, 'plugin_title_scaled')
                set(th, 'FontSize', get(th, 'FontSize') * titlescale);
                setappdata(th, 'plugin_title_scaled', true);
            end

        end

    end

end

end % function
