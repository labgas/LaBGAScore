function outFiles = save_all_open_figures_smart(outDir, prefix, formats, closeFigs)
% save_all_open_figures_smart
%
% Save all open MATLAB figures with smart filenames.
%
% Priority for naming:
%   1) figure Name
%   2) axes title
%   3) fallback: fig01, fig02, ...
%
% INPUTS
%   outDir     - output directory (default: pwd)
%   prefix     - filename prefix (default: 'Figure')
%   formats    - cell array, e.g. {'fig','svg'} (default: {'fig','svg'})
%   closeFigs  - logical (default: false)
%
% OUTPUT
%   outFiles   - struct with saved file paths
%
% EXAMPLE
%   save_all_open_figures_smart('results','PLSDA',{'fig','svg'},false)

if nargin < 1 || isempty(outDir)
    outDir = pwd;
end
if nargin < 2 || isempty(prefix)
    prefix = 'Figure';
end
if nargin < 3 || isempty(formats)
    formats = {'fig','svg'};
end
if nargin < 4
    closeFigs = false;
end

if ischar(formats) || isstring(formats)
    formats = cellstr(formats);
end

if ~exist(outDir,'dir')
    mkdir(outDir);
end

figs = findall(0,'Type','figure');

if isempty(figs)
    warning('No open figures to save.');
    outFiles = struct([]);
    return
end

% Sort by figure number
[~,idx] = sort([figs.Number]);
figs = figs(idx);

outFiles = struct('figureNumber',{},'name',{},'files',{});
usedNames = containers.Map('KeyType','char','ValueType','double');

for i = 1:numel(figs)

    h = figs(i);

    % ---------- Get smart name ----------
    name = get(h,'Name');

    if isempty(name)
        name = get_axes_title(h);
    end

    if isempty(name)
        name = sprintf('fig%02d',i);
    end

    name = sanitize_filename(name);
    prefixSafe = sanitize_filename(prefix);

    % Avoid duplicates
    if isKey(usedNames,name)
        usedNames(name) = usedNames(name) + 1;
        name = sprintf('%s_%02d', name, usedNames(name));
    else
        usedNames(name) = 1;
    end

    base = fullfile(outDir, sprintf('%s_%s', prefixSafe, name));

    saved = {};

    % ---------- Save formats ----------
    for f = 1:numel(formats)

        fmt = lower(formats{f});

        switch fmt
            case 'fig'
                file = [base '.fig'];
                savefig(h,file,'compact');

            case 'svg'
                file = [base '.svg'];
                try
                    exportgraphics(h,file,'ContentType','vector');
                catch
                    print(h,file,'-dsvg');
                end

            case 'png'
                file = [base '.png'];
                exportgraphics(h,file,'Resolution',300);

            case 'pdf'
                file = [base '.pdf'];
                exportgraphics(h,file,'ContentType','vector');

            otherwise
                warning('Unsupported format: %s',fmt);
                continue
        end

        saved{end+1} = file;
    end

    outFiles(i).figureNumber = h.Number;
    outFiles(i).name = name;
    outFiles(i).files = saved;

    if closeFigs
        close(h);
    end
end

fprintf('Saved %d figure(s) to %s\n', numel(figs), outDir);

end


% =========================
% Helper: get axes title
% =========================
function t = get_axes_title(hFig)

t = '';

ax = findall(hFig,'Type','axes');

for k = 1:numel(ax)
    try
        titleStr = get(get(ax(k),'Title'),'String');
    catch
        titleStr = '';
    end

    if iscell(titleStr)
        titleStr = strjoin(titleStr,' ');
    end

    if isstring(titleStr)
        titleStr = char(titleStr);
    end

    titleStr = strtrim(titleStr);

    if ~isempty(titleStr)
        t = titleStr;
        return
    end
end
end


% =========================
% Helper: sanitize filename
% =========================
function s = sanitize_filename(s)

if isstring(s)
    s = char(s);
end

s = strtrim(s);

% Replace spaces
s = regexprep(s,'\s+','_');

% Remove problematic chars
s = regexprep(s,'[^\w\-\.]','_');

% Collapse underscores
s = regexprep(s,'_+','_');

% Trim
s = regexprep(s,'^_+|_+$','');

if isempty(s)
    s = 'figure';
end

end

