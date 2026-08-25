function setParforStream(seed, stageOffset, iter)
% setParforStream  Make one parfor iteration's random draws reproducible.
%
%   setParforStream(seed, stageOffset, iter)
%
%   Installs a Threefry random stream positioned at substream
%   (stageOffset + iter) as the global stream for the current iteration.
%   The draws an iteration makes then depend only on (seed, stageOffset, iter)
%   -- not on which worker ran it, nor on how many workers the pool has, nor on
%   the order iterations happen to complete in.
%
%   Inputs
%     seed         base RNG seed (opts.seed)
%     stageOffset  large constant separating pipeline stages, so that the outer
%                  CV, permutation, bootstrap and learning-curve loops never
%                  draw from the same substreams
%     iter         the parfor loop index
%
%   Why the global stream is set rather than a stream object passed around
%     randperm and randsample accept an explicit stream argument, but
%     cvpartition does not, and every stage here builds partitions. Setting the
%     global stream is therefore unavoidable if fold assignment is to be
%     reproducible.
%
%   What this fixes
%     All four pipelines seeded once on the client (rng(1) or rng(opts.seed)),
%     but their permutation, bootstrap and learning-curve loops are parfor loops
%     whose randperm/randsample/cvpartition calls execute on workers, whose
%     streams the client seed does not reach. Those stages were therefore not
%     reproducible from run to run. For ENet the outer repeat loop is itself a
%     parfor containing cvpartition, so even its headline AUC was not
%     reproducible.
%
%   Stage offsets in use across the family:
%     1e6  outer repeated CV
%     2e6  permutation
%     3e6  bootstrap
%     4e6  learning curve
%
%   See also RandStream, parfor, cvpartition.

s = RandStream('Threefry','Seed',seed);
s.Substream = stageOffset + iter;
RandStream.setGlobalStream(s);

end
