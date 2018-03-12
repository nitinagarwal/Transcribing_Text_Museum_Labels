function [newpoints,error,coeff] = linearFit(points)
%fit a straight line

% y = Ax+B

x=points(:,1);
y=points(:,2);

A=[x(:).^1 ones(size(x))];
coeff = (A' * A)\(A' * y);
newpoints = [x A*coeff];

error= sum((newpoints(:,1)-points(:,1)).^2 + (newpoints(:,2)-points(:,2)).^2);
