function [mag_histogram] = get_histogram2(mag, binning, nBin)
% instead of computing the mag and ori histogram separate
% we compute the average mag : total mag/ total number of points with that mag

mag_histogram = zeros(nBin, 1);  % histogram of both the mag and ori
ori_histogram = zeros(nBin, 1);

for i = 1:nBin
    mag_histogram(i) = sum(mag(binning == i));
    ori_histogram(i) = numel(find(binning == i));    
end;


mag_histogram = mag_histogram./ori_histogram;

mag_histogram(isnan(mag_histogram))=0;  % for some orientation number could be 0. 


