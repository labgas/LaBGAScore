%% LaBGAScore_smart_parallel_pool_setup.m
%
%
% *USAGE*
%
% This script sets up (or resizes) a Matlab parallel pool sized to a
% capped fraction of available workers, so that parallel jobs don't
% claim every core on a shared machine. Uses 66% of currently available
% workers, always leaving at least 1 worker free overall.
%
% -------------------------------------------------------------------------
%
% modified by: Lukas Van Oudenhove
%
% date:   April, 2026
%
% -------------------------------------------------------------------------
%
%
%% SET UP PARALLEL POOL
% -------------------------------------------------------------------------

pc = parcluster;
maxWorkers = pc.NumWorkers;

pool = gcp('nocreate');

if isempty(pool)
    usedWorkers = 0;
else
    usedWorkers = pool.NumWorkers;
end

availableWorkers = maxWorkers - usedWorkers;

% Use 66% of available workers
% 0.75, raised from 0.66 on 2026-09-25. The cap exists so concurrent jobs do
% not oversubscribe the box; it works only if the PROFILE NumWorkers reflects
% PHYSICAL cores. On this machine (AMD EPYC 7282, 2 sockets x 16 = 32 physical,
% 64 logical via SMT2) the profile is 32, so 0.66 gave 21 and 0.75 gives 24 -
% both comfortably inside the physical count.
%
% DO NOT set the profile from /proc's logical count, and do not call saveProfile
% to raise it. Both were done here on 2026-09-25 (profile pushed to 40) and a
% permutation job then ran ~7x slower than predicted: every worker pegged at
% 99%, but 8 of them on hyperthreads and all of them thrashing a 64 MiB
% per-socket L3 with a 175 MB feature matrix. The cap is the safeguard; sizing
% the profile wrongly defeats it.
nWorkers = floor(0.75 * availableWorkers);

% Cap: always leave at least 1 worker free overall
nWorkers = min(nWorkers, maxWorkers - 1);

% Ensure at least 1 worker
nWorkers = max(1, nWorkers);

% Start or resize pool if needed
if isempty(pool) || pool.NumWorkers ~= nWorkers
    if ~isempty(pool)
        delete(pool);
    end
    parpool(pc, nWorkers);
end

fprintf('Using %d/%d workers (%.0f%% of available, capped)\n', ...
    nWorkers, maxWorkers, 100 * nWorkers / maxWorkers);