function report = LaBGAScore_stats_rederive_storey_q(modeldirs, varargin)
%% LaBGAScore_stats_rederive_storey_q.m
%
%
% *USAGE*
%
% Re-derives Storey q-values in ALREADY SAVED second-level results, without
% re-running any analysis, and reports what changes.
%
% This exists because q = pi0 * q_BH (floored at p) is a scalar transform of the
% p-values, and the p-values are stored in the results tables next to the
% q-values. Everything expensive - model fitting, bootstrapping, permutation,
% image loading - produced those p-values; Storey is applied afterwards. So when
% the pi0 estimator changes, the q-values can be corrected by reading the saved
% tables. There is no need to re-run prep_3a, c2a or the decoding scripts.
%
% The function executes the following steps, per model directory:
%
% 1. find every <modeldir>/results/*stats*.mat
%
% 2. search each file for tables carrying BOTH a 'p' and a 'q_Storey' column,
%   including tables nested inside cells and struct fields (roi_glm_stats{c},
%   parcelwise and neurotransmitter stats are all shaped that way)
%
% 3. recompute pi0 and q with the CURRENT LaBGAScore_Storey_FDR, and compare
%   against what is stored
%
% 4. report, per table: n, max(p), old and new pi0, the largest change in q, how
%   many results sat below alpha before and after, how many crossed alpha, and
%   whether the stored q was simply the raw p-values
%   NOTE: q == p happens when pi0 is small enough that pi0*q_BH falls below p
%       and the mandatory q >= p floor takes over. The column is then an
%       UNCORRECTED p-value under an FDR name - worth knowing independently of
%       any pi0 change
%
% 5. with 'write' true, save corrected copies of the tables and the report
%   NOTE: originals are NEVER modified. Corrected tables go to a new file with
%       a '_storeyfix' tag beside the original, so the old numbers remain
%       available for comparison and nothing already published is overwritten
%
% Report-only is the default, deliberately: run it, read the diff, then decide
% per model whether anything needs re-publishing.
%
%
% *OPTIONS*
%
% * alpha           significance threshold used only for the crossing counts, default 0.05
%
% * write           default false; true also writes corrected tables and the report to disk
%
% * outdir          where the report goes, default <modeldir>/results/notes
%
% * method          method passed to LaBGAScore_Storey_FDR, default 'sas'
%
% * verbose         default true; print the per-table summary as it goes
%
%
% *DEPENDENCIES*
%
% 1. LaBGAScore Github repo on Matlab path, with subfolders
%   https://github.com/labgas/LaBGAScore
%   in particular stats_tools/functions/LaBGAScore_Storey_FDR.m
% 2. No CanlabCore or SPM dependency: this reads .mat files and tables only
%
%
% *NOTES*
%
% INPUTS: a second-level model directory, or a cell array of them, e.g.
%   {'/data/proj_cfs/secondlevel/model_2c_IOM', ...
%    '/data/proj_discoverie/secondlevel/model_2h_casecontrol_zcond_combat'}
%
% OUTPUT: a table, one row per recomputed stats table, also written as
% storey_rederive_report.tsv under outdir when 'write' is true
%
% WHAT THIS DOES NOT COVER. Storey is also applied to large-n voxelwise and
% parcelwise p-value vectors inside prep_3a and the decoding scripts. Those
% vectors are not always saved as tables, and at large n the upper tail of the
% lambda curve is populated, which is exactly where pi0 estimators agree - so
% the expected change there is small. Confirm with this report first; only
% re-run prep_3a if it shows movement.
%
% -------------------------------------------------------------------------
%
% by: Lukas Van Oudenhove  |  KU Leuven, October 2026
%
% -------------------------------------------------------------------------
%
% LaBGAScore_stats_rederive_storey_q.m         v1.0
%
% last modified: 2026/10/01
%
%

%% PARSE OPTIONS
% -------------------------------------------------------------------------

ip = inputParser;
ip.addParameter('alpha',   0.05,  @isnumeric);
ip.addParameter('write',   false, @(x) islogical(x) || isnumeric(x));
ip.addParameter('outdir',  '',    @ischar);
ip.addParameter('method',  'sas', @(x) ischar(x) || isstring(x));
ip.addParameter('verbose', true,  @(x) islogical(x) || isnumeric(x));
ip.parse(varargin{:});
alpha   = ip.Results.alpha;
dowrite = logical(ip.Results.write);
outdir  = ip.Results.outdir;
method  = char(ip.Results.method);
verbose = logical(ip.Results.verbose);

if ischar(modeldirs) || isstring(modeldirs)
    modeldirs = {char(modeldirs)};
end

if isempty(which('LaBGAScore_Storey_FDR'))
    error('LaBGAScore_stats_rederive_storey_q:missingDep', ...
        '\nLaBGAScore_Storey_FDR not found on the Matlab path, please add the LaBGAScore repo WITH subfolders');
end

rows = {};

%% LOOP OVER MODELS AND RESULTS FILES
% -------------------------------------------------------------------------

for d = 1:numel(modeldirs)

    modeldir = modeldirs{d};
    [~, modelname] = fileparts(modeldir);
    resdir = fullfile(modeldir, 'results');

        if ~isfolder(resdir)
            warning('LaBGAScore_stats_rederive_storey_q:noResults', ...
                '\nno results dir in %s, skipping', modeldir);
            continue
        end

    files = dir(fullfile(resdir, '*stats*.mat'));

        if verbose
            fprintf('\n%s\n  %d results file(s) with "stats" in the name\n', modelname, numel(files));
        end

    for f = 1:numel(files)

        fpath = fullfile(files(f).folder, files(f).name);
        S = load(fpath);
        changed_any = false;
        Sout = S;

        fn = fieldnames(S);
        for v = 1:numel(fn)

            X = S.(fn{v});
            % tables arrive bare, or wrapped one level deep in a cell - both occur
            wascell = iscell(X);
                if ~wascell
                    X = {X};
                end

            for c = 1:numel(X)

                T = X{c};
                    if ~istable(T), continue, end
                cols = T.Properties.VariableNames;
                    if ~all(ismember({'p','q_Storey'}, cols)), continue, end

                p  = T.p(:);
                ok = ~isnan(p);
                    if sum(ok) < 4
                        % below 4 p-values no pi0 estimator means anything
                        continue
                    end

                [q_new, pi0_new, info] = LaBGAScore_Storey_FDR(p(ok), ...
                    'method', method, 'verbose', false);

                q_old   = T.q_Storey(ok);
                pi0_old = NaN;
                    if ismember('pi0', cols), pi0_old = T.pi0(find(ok,1)); end

                % stored q identical to p means the q >= p floor took over, i.e.
                % the column is an uncorrected p-value under an FDR name
                q_was_p = all(abs(q_old - p(ok)) < 1e-10);

                rows(end+1,:) = { modelname, files(f).name, ...
                    sprintf('%s{%d}', fn{v}, c), sum(ok), max(p(ok)), ...
                    pi0_old, pi0_new, max(abs(q_new - q_old)), ...
                    sum(q_old < alpha), sum(q_new < alpha), ...
                    sum((q_old < alpha) ~= (q_new < alpha)), ...
                    q_was_p, info.reliable, strjoin(info.reasons, '; ') }; %#ok<AGROW>

                    if verbose
                        fprintf(['    %-30s %-18s n=%-4d pi0 %.4f -> %.4f  ' ...
                                 'dq<=%.4f  <%.2f: %d -> %d  flips %d%s\n'], ...
                            files(f).name(1:min(30,end)), sprintf('%s{%d}', fn{v}, c), ...
                            sum(ok), pi0_old, pi0_new, max(abs(q_new - q_old)), ...
                            alpha, sum(q_old < alpha), sum(q_new < alpha), ...
                            sum((q_old < alpha) ~= (q_new < alpha)), ...
                            repmat('   [stored q WAS raw p]', 1, q_was_p));
                    end

                    if dowrite
                        Tnew = T;
                        Tnew.q_Storey(ok) = q_new;
                            if ismember('pi0', cols),             Tnew.pi0(ok) = pi0_new; end
                            if ismember('storey_reliable', cols), Tnew.storey_reliable(ok) = info.reliable; end
                        X{c} = Tnew;
                        changed_any = true;
                    end
            end

                if dowrite && changed_any
                    if wascell, Sout.(fn{v}) = X; else, Sout.(fn{v}) = X{1}; end
                end
        end

            if dowrite && changed_any
                [~, base, ext] = fileparts(files(f).name);
                newpath = fullfile(files(f).folder, [base '_storeyfix' ext]);
                save(newpath, '-struct', 'Sout', '-v7.3');
                    if verbose
                        fprintf('    -> wrote %s\n', [base '_storeyfix' ext]);
                    end
            end
    end
end

%% BUILD AND OPTIONALLY WRITE THE REPORT
% -------------------------------------------------------------------------

varnames = {'model','file','table','n','max_p','pi0_old','pi0_new','dq_max', ...
            'n_sig_old','n_sig_new','n_cross_alpha','stored_q_was_raw_p', ...
            'reliable_new','reasons_new'};

if isempty(rows)
    report = cell2table(cell(0, numel(varnames)), 'VariableNames', varnames);
    warning('LaBGAScore_stats_rederive_storey_q:nothingFound', ...
        '\nno saved table with both a p and a q_Storey column was found - nothing to re-derive');
    return
end

report = cell2table(rows, 'VariableNames', varnames);

if verbose
    fprintf('\n%d table(s) re-derived\n', height(report));
    fprintf('  tables whose q changed by > 0.001 : %d\n', sum(report.dq_max > 0.001));
    fprintf('  tables crossing alpha = %.2f       : %d\n', alpha, sum(report.n_cross_alpha > 0));
    fprintf('  tables whose stored q was raw p    : %d\n', sum(report.stored_q_was_raw_p));
    fprintf('  tables with an unreliable new pi0  : %d\n', sum(~report.reliable_new));
    if ~dowrite
        fprintf('\n  REPORT ONLY - nothing written. Re-run with ''write'', true to save\n');
        fprintf('  corrected tables (as *_storeyfix.mat, originals untouched).\n');
    end
end

if dowrite
    if isempty(outdir)
        outdir = fullfile(modeldirs{1}, 'results', 'notes');
    end
    if ~isfolder(outdir), mkdir(outdir); end
    tsv = fullfile(outdir, 'storey_rederive_report.tsv');
    writetable(report, tsv, 'FileType', 'text', 'Delimiter', 'tab');
    if verbose, fprintf('\n  report written to %s\n', tsv); end
end

end
