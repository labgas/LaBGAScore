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
nWorkers = floor(0.66 * availableWorkers);

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