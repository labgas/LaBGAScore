function p = pvals_from_ranks(x)
% Uncorrected, voxel-wise empirical p-values

x = x(:);
n = numel(x);

[~,ord] = sort(x, 'descend');
r = zeros(n,1);
r(ord) = 1:n;

p = r ./ (n + 1);

end

