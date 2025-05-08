% Parameters
base_corr = 0.85;
decay_rate = 1; % This doesn't matter when decay rate is 1 — simplifies to abs(i-j)
n = 5; % Matrix size

% Initialize LEAR matrix
LEAR = eye(n);

% Fill off-diagonal entries using exponential decay
for i = 1:n
    for j = 1:n
        if i ~= j
            LEAR(i,j) = base_corr ^ abs(i-j);
        end
    end
end

% Display the result
disp('LEAR Matrix:');
disp(LEAR);