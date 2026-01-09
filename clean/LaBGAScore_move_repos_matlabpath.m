function LaBGAScore_move_repos_matlabpath(repoA, repoB)

% moves repoA below repoB in matlab path and saves new path - saves you
% tons of clicking in pathtool ;)

p = strsplit(path, pathsep);

idxA = find(contains(p, repoA));

block = p(idxA);
p(idxA) = [];

idxB = find(contains(p, repoB), 1, 'last');

p = [p(1:idxB) block p(idxB+1:end)];

path(strjoin(p, pathsep))
savepath
path

end