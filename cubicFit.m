function [newpoints,error,coeff] = cubicFit(points)
%fit a cubic polynomial
% y = Ax^3+Bx^2+Cx+D

x=points(:,1);
y=points(:,2);

A=[x(:).^3 x(:).^2 x(:).^1 ones(size(x))];
coeff = (A' * A)\(A' * y);
newpoints = [x A*coeff];

error= sum(sqrt((newpoints(:,1)-points(:,1)).^2 + (newpoints(:,2)-points(:,2)).^2));
