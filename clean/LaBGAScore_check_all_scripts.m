function issues = LaBGAScore_check_all_scripts(rootdir)
%
%
% *USAGE*
%
% Runs MATLAB's built-in Code Analyzer (checkcode) on every .m file found
% recursively under rootdir, and prints a report to the command window.
%
% Syntax errors (checkcode id SYNER) are reported separately and
% prominently. Tested against real bugs found and fixed in this repo:
% checkcode reliably catches genuine parse errors (a file that cannot run
% past that point, e.g. a stray "sort{}"), but does NOT catch undefined
% variables used at runtime, calls to functions that don't exist or
% aren't on the path, or logic bugs (e.g. wrong array indexing) - those
% were only found by actually reading the code against its documented
% behavior. Treat this as a fast baseline safety net, not a substitute
% for review.
%
% All other Code Analyzer findings (style/performance suggestions,
% unreachable code, etc.) are summarized per file rather than printed in
% full, since most of them are not actionable bugs; use checkcode
% directly on a given file for full detail.
%
%
% *OPTIONS*
%
% * rootdir     directory tree to check, default: this repo's own root,
%               determined from this function's own location (i.e. one
%               level up from the clean/ folder it lives in)
%
%
% -------------------------------------------------------------------------
%
% modified by: Lukas Van Oudenhove & Claude Sonnet 5
%
% date:   KU Leuven, August, 2026
%
% -------------------------------------------------------------------------
%
% LaBGAScore_check_all_scripts.m         v1.0
%
% last modified: 2026/08/20
%
%
%% SET DEFAULT ROOTDIR
% -------------------------------------------------------------------------

if nargin < 1 || isempty(rootdir)
    rootdir = fileparts(fileparts(mfilename('fullpath'))); % clean/ -> repo root
end


%% FIND ALL .m FILES AND RUN CODE ANALYZER
% -------------------------------------------------------------------------

filelist = dir(fullfile(rootdir,'**','*.m'));
filelist = filelist(~[filelist.isdir]);

file = strings(0,1);
line = zeros(0,1);
id = strings(0,1);
message = strings(0,1);

for f = 1:numel(filelist)

    fullpath = fullfile(filelist(f).folder, filelist(f).name);
    relpath = erase(fullpath, [rootdir filesep]);

    info = checkcode(fullpath,'-id');

    for k = 1:numel(info)
        file(end+1,1) = string(relpath); %#ok<AGROW>
        line(end+1,1) = info(k).line;
        id(end+1,1) = string(info(k).id);
        message(end+1,1) = string(info(k).message);
    end

end

issues = table(file,line,id,message,'VariableNames',{'File','Line','Id','Message'});
issues.IsSyntaxError = issues.Id == "SYNER";


%% PRINT REPORT
% -------------------------------------------------------------------------

fprintf('\n%d .m files checked, %d total Code Analyzer findings\n', numel(filelist), height(issues));

synerrors = issues(issues.IsSyntaxError,:);

if isempty(synerrors)
    fprintf('No syntax errors found.\n\n');
else
    fprintf('\n*** %d SYNTAX ERROR(S) - these files will not run past the flagged line ***\n', height(synerrors));
    for k = 1:height(synerrors)
        fprintf('  %s (line %d): %s\n', synerrors.File(k), synerrors.Line(k), synerrors.Message(k));
    end
    fprintf('\n');
end

other = issues(~issues.IsSyntaxError,:);

if ~isempty(other)
    [countsPerFile, filenames] = groupcounts(other.File);
    [countsPerFile, sortidx] = sort(countsPerFile,'descend');
    filenames = filenames(sortidx);
    fprintf('Other Code Analyzer findings per file (style/performance suggestions - not verified to catch undefined variables, wrong function calls, or logic bugs; review at your discretion):\n');
    for k = 1:numel(filenames)
        fprintf('  %-70s %d\n', filenames(k), countsPerFile(k));
    end
    fprintf('\nRun checkcode(''<file>'',''-id'') directly on any file above for full detail.\n\n');
end

end
