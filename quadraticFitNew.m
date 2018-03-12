function [newpoints,error,coeff] = quadraticFitNew(points)
%fit a cubic polynomial

% y = Ax^2+Bx+C

x=points(:,1);
y=points(:,2);

A=[x(:).^2 x(:).^1 ones(size(x))];
coeff = (A' * A)\(A' * y);
newpoints = [x A*coeff];

error= sum((newpoints(:,1)-points(:,1)).^2 + (newpoints(:,2)-points(:,2)).^2);

% Method II. 
% % kernel regression...
% % scale x down for numerical stability...
% sigma = 60;
% lambda = 5;
% 
% % x = traj(:, 1);
% % y = traj(:, 2);
% 
% xscale = max(abs(x));
% x = x / xscale;
% sigma = sigma / xscale;
% 
% 
% % assume f(x) = ax^3 + bx^2 + cx + d + \sum_i c_i * rbf(x_i, x)
% n = length(x);
% dists = abs(repmat(x(:), 1, n) - repmat(x(:)', n, 1));
% rbf = exp(-dists.^2 / 2 / sigma.^2);
% 
% % minimize the following function...
% % min \sum_j (y_j - f(x_j))^2 + lambda * ||c||^2 
% % or min (A[a, b, c, d, ci] - y).^2 + lambda c'*c
% 
% A = [x(:).^3, x(:).^2, x(:), ones(n, 1), rbf];
% coeffs = (A'*A + diag([0, 0, 0, 0, repmat(lambda, 1, n)])) \ (A' * y);
% 
% newpoints = [points(:,1), A * coeffs];
% 
% coeffs(1:4) = coeffs(1:4) ./ [xscale.^3; xscale.^2; xscale; 1];


end