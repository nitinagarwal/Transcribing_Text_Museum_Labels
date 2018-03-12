function [magSqr, binning,ori] = im_bin(gx, gy, nBin)
% NOTE: where 0 is starting from (from the left and not from the right)

% For normal gx and gy: 
%  (1,1)
%   --------------------------------------------------------------
%   | Image
%   |
%   |                        nBin/4 + 1
%   |                              |
%   |                              |
%   |                              |
%   |                1 <---------- O --------->  nBin/2+1
%   |                              |              
%   |                              |           
%   |                              |       ...
%   |                       3*nBin/4 + 1

[m, n] = size(gx);

magSqr = sqrt(gx.^2 + gy.^2);

% range: (-180, 180)
ori = atan2d(gy, gx);
%range : (0,180)
ori((ori<0)) = ori((ori<0))+180;

edges = linspace(0, 180, nBin + 1);
edgesBoundary = (edges(1:end-1) + edges(2:end)) / 2;

binning = ones(m, n);
for i = 1:nBin-1
    binning = binning + (ori > edgesBoundary(i));
end; 
binning(ori > edgesBoundary(nBin)) = 1;


