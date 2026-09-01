function R = LaBGAScore_check_display(varargin)
%
%
% *USAGE*
%
% Checks whether the display of the current session is big enough to
% produce good figures in a published HTML report, and says what to change
% if it is not.
%
% Why this exists: publish() captures figures from the screen, so the size
% and DPI of your X2go session decide how the figures in your report come
% out. Two things follow, and neither is obvious:
%
% # A figure larger than the session window cannot be captured at the size
%   it was asked for. plugin_set_figure_size scales it down to fit,
%   preserving the aspect ratio, so nothing is distorted - but the figure
%   is smaller than intended, and differs from a colleague's.
% # Figure size is specified in INCHES, not pixels, because MATLAB font
%   sizes are in points (1 pt = 1/72 inch). Inches available = pixels / DPI,
%   so RAISING the session DPI REDUCES the figure size you can achieve at a
%   given window size. A high DPI makes menus comfortable and figures
%   small; they trade against each other.
%
%
% *OPTIONS*
%
% * 'target'    [width height] inches you want figures to reach. Default
%               [12 7.5], which is plugin_set_figure_size's default and is
%               reachable on every screen in the lab at 96 DPI. Keep the two
%               in step if you change either.
% * 'margin'    [horizontal vertical] screen fraction reserved for window
%               decorations, default [0.02 0.06] (matches the plugin)
% * 'print'     logical, default true
%
%
% *OUTPUT*
%
% * R           struct with the measured display, the largest figure it can
%               produce, whether the target fits, the single DPI value to
%               set (R.recommended_dpi) and what else to change
%
%
% *NOTES*
%
% X2go's "custom" width/height is only the size the session STARTS at. If
% you resize the X2go window afterwards, the X display is resized with it
% and the new size is what counts. This function reports what is true now,
% which is what publish() will actually capture.
%
%
% *SEE ALSO*
%
% plugin_set_figure_size (CANlab_help_examples), LaBGAScore_prov_publish,
% and the "Display settings for publishing figures" section of
% LaBGAS_fMRI_analysis_workflow.md
%
% -------------------------------------------------------------------------
%
% modified by: Lukas Van Oudenhove & Claude Opus 5
%
% date:   KU Leuven, August, 2026
%
% -------------------------------------------------------------------------
%
% LaBGAScore_check_display.m         v1.0
%
% last modified: 2026/08/31
%
%
%% PARSE OPTIONS
% -------------------------------------------------------------------------

p = inputParser;
p.addParameter('target', [12 7.5], @(x) isnumeric(x) && numel(x) == 2 && all(x > 0));
p.addParameter('margin', [0.02 0.06], @(x) isnumeric(x) && numel(x) == 2);
p.addParameter('print', true, @islogical);
p.parse(varargin{:});
opt = p.Results;

target = opt.target;
margin = opt.margin;


%% MEASURE THE DISPLAY
% -------------------------------------------------------------------------

ss = get(0, 'ScreenSize');
dpi = get(0, 'ScreenPixelsPerInch');

R.screen_px = [ss(3) ss(4)];
R.dpi = dpi;
R.usable_px = floor([ss(3)*(1-margin(1)), ss(4)*(1-margin(2))]);
R.usable_in = R.usable_px / dpi;
R.target_in = target;
R.target_px = ceil(target * dpi);

scale = min([1, R.usable_in(1)/target(1), R.usable_in(2)/target(2)]);
R.achievable_in = target * scale;
R.achievable_px = round(R.achievable_in * dpi);
R.fits = scale >= 1;

% at this window size, the highest DPI that would let the target fit
R.max_dpi_for_target = floor(min(R.usable_px ./ target));

% A single value to set, rather than a ceiling. Provided the target fits,
% a HIGHER DPI is better: the figure stays the same physical size, so the
% font-to-canvas ratio is unchanged, but it is captured at more pixels and
% is therefore sharper in the report - and on-screen text is more
% comfortable. So recommend the highest value people actually use that
% still leaves a little headroom for panels and window decorations.
standard_dpi = [96 110 120 133 144];
candidates = standard_dpi(standard_dpi <= R.max_dpi_for_target * 0.97);
if isempty(candidates)
    R.recommended_dpi = 96;                 % nothing fits; 96 is the floor
else
    R.recommended_dpi = max(candidates);
end
R.recommended_figure_px = round(target * R.recommended_dpi);

% at this DPI, the window size the target would need
R.min_window_px = ceil([target(1)*dpi/(1-margin(1)), target(2)*dpi/(1-margin(2))]);


%% PRINT A REPORT AND A RECOMMENDATION
% -------------------------------------------------------------------------

if ~opt.print
    return
end

dashes = repmat('-', 1, 74);

fprintf('\n%s\nDISPLAY CHECK FOR PUBLISHED FIGURES\n%s\n', dashes, dashes);

fprintf('session window   %d x %d px at %g DPI\n', R.screen_px(1), R.screen_px(2), dpi);
fprintf('usable for a figure  %d x %d px  =  %.1f x %.1f inches\n', ...
    R.usable_px(1), R.usable_px(2), R.usable_in(1), R.usable_in(2));
fprintf('target figure    %g x %g inches  (needs %d x %d px)\n\n', ...
    target(1), target(2), R.target_px(1), R.target_px(2));

if R.fits
    fprintf('GOOD: this session can produce the full %g x %g inch figure.\n', target);
    fprintf('Reports made here will match anyone else whose session also fits it.\n');

    if R.recommended_dpi > dpi
        fprintf('\nOptional: raising the DPI to %d would still fit, and would capture\n', ...
            R.recommended_dpi);
        fprintf('the same %g x %g inch figure at %d x %d px instead of %d x %d - sharper in\n', ...
            target(1), target(2), R.recommended_figure_px(1), R.recommended_figure_px(2), ...
            round(target(1)*dpi), round(target(2)*dpi));
        fprintf('the report, and larger on-screen text. Figures stay the same physical\n');
        fprintf('size, so they remain comparable with colleagues.\n');
    elseif R.recommended_dpi < dpi
        fprintf('\nYour DPI (%d) is above the %d this window can comfortably support,\n', ...
            dpi, R.recommended_dpi);
        fprintf('but still fits. Nothing to change unless figures start being scaled.\n');
    else
        fprintf('\nYour DPI (%d) is the recommended setting for this window size.\n', dpi);
    end
else
    fprintf('TOO SMALL: figures will be scaled down to %.1f x %.1f inches\n', ...
        R.achievable_in(1), R.achievable_in(2));
    fprintf('(%d x %d px). Aspect ratio is preserved, so nothing is distorted,\n', ...
        R.achievable_px(1), R.achievable_px(2));
    fprintf('but text will look relatively larger than in a colleague''s report.\n\n');

    fprintf('EITHER set the session DPI to %d, keeping this window size\n', ...
        R.recommended_dpi);
    fprintf('       (X2go: Session preferences > Input/Output > "Set display DPI")\n');
    fprintf('OR     enlarge the X2go window to at least %d x %d px at %g DPI\n', ...
        R.min_window_px(1), R.min_window_px(2), dpi);
    fprintf('       (just drag the window; X2go resizes the session with it)\n');

    if R.max_dpi_for_target < 72
        fprintf('\nNote: %d DPI is impractically low. This window is too small for a\n', ...
            R.max_dpi_for_target);
        fprintf('%g x %g inch figure at any sensible DPI - enlarge the window instead,\n', target);
        fprintf('or agree a smaller lab-wide target (see the workflow documentation).\n');
    end
end

fprintf('\n%s\nWHAT YOUR WINDOW SIZE ALLOWS, BY DPI\n%s\n', dashes, dashes);
fprintf('  %-8s %-18s %s\n', 'DPI', 'LARGEST FIGURE', 'FITS TARGET?');
for d = [72 96 110 120 133 150]
    inches = R.usable_px / d;
    ok = all(inches >= target);
    if ok, verdict = 'yes'; else, verdict = sprintf('no (max %.1f x %.1f in)', inches); end
    marker = ' ';
    if d == dpi, marker = '*'; end
    fprintf('%s %-8d %-18s %s\n', marker, d, sprintf('%.1f x %.1f in', inches), verdict);
end
fprintf('  (* = your current setting)\n');

fprintf('\n%s\nREMINDER\n%s\n', dashes, dashes);
fprintf('X2go''s "custom" width/height is only the size the session STARTS at.\n');
fprintf('Resizing the X2go window resizes the session, and the new size is what\n');
fprintf('publish() captures. The numbers above are what is true right now.\n');
fprintf('%s\n\n', dashes);

end
