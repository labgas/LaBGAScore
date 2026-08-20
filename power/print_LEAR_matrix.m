%% print_LEAR_matrix.m
%
%
% *USAGE*
%
% Small demo script that builds and displays a Linear Exponential
% AutoRegressive (LEAR)-style correlation matrix, in which off-diagonal
% entries decay exponentially with distance from the diagonal
% (base_corr ^ abs(i-j))
%
%
% *OPTIONS*
%
% * base_corr        base correlation used for the exponential decay, default 0.85
% * n                 size of the (square) LEAR matrix, default 5
%
%
% *NOTES*
%
% * decay_rate is defined in the code but is currently a no-op - the matrix
%   formula (base_corr ^ abs(i-j)) never references it, so changing its
%   value has no effect
%
%
% -------------------------------------------------------------------------
%
% modified by: Lukas Van Oudenhove
%
% date:   May, 2025
%
% -------------------------------------------------------------------------
%
% print_LEAR_matrix.m
%
% last modified: 2026/08/20
%
%
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