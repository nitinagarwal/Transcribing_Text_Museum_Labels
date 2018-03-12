function newpoints = linearFit_evaluate(xvalues,coeffs)
% coeff is 2x1

if(size(xvalues,2)~=1)
    xvalues=xvalues';
end

if(size(coeffs,2)~=1)
    coeffs=coeffs';
end

x=xvalues;

A=[x(:).^1 ones(size(x))];
coeff = coeffs;

newpoints = [xvalues A*coeff];


end