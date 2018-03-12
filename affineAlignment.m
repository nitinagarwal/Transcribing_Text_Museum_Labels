function [I1, newimage,tform] = affineAlignment(I1,I2)
% warping I2 

if(ndims(I1)==3)
    im1=rgb2gray(I1);
end

if(ndims(I2)==3)
    im2=rgb2gray(I2);
end

%-------feature detection & extraction using SURF
im1_pts = detectSURFFeatures(im1);
im2_pts = detectSURFFeatures(im2);

% im1_pts = detectHarrisFeatures(im1,'MinQuality',.001);
% im2_pts = detectHarrisFeatures(im2,'MinQuality',.001);

[featuresSURF_im1,validptsSURF_im1] = extractFeatures(im1,im1_pts,'Method','SURF','SURFSize',128);
[featuresSURF_im2,validptsSURF_im2] = extractFeatures(im2,im2_pts,'Method','SURF','SURFSize',128);

[im1SIFT_pts,featuresSIFT_im1] = vl_sift(single(im1));
[im2SIFT_pts,featuresSIFT_im2] = vl_sift(single(im2));

%-------matching-------
[indexPairs,~] = matchFeatures(featuresSURF_im1,featuresSURF_im2,'MatchThreshold',20,'MaxRatio',.8);
matched_im1  = validptsSURF_im1(indexPairs(:,1));
matched_im2 = validptsSURF_im2(indexPairs(:,2));

% outlier removal
[~, inlier_im2,inlier_im1] = estimateGeometricTransform(matched_im2,matched_im1,'affine','MaxDistance',20);

% showCorresponces(I1,I2,matched_im1.Location(),matched_im2.Location(),true)
% showCorresponcesVertical(im1,im2,matched_im1.Location,matched_im2.Location,false)
% figure,showMatchedFeatures(im1,im2,matched_im1,matched_im2,'montage')


[matches, scores] = vl_ubcmatch(featuresSIFT_im1, featuresSIFT_im2,5) ;
matched_im1 = (im1SIFT_pts(1:2,matches(1,:)))';
matched_im2 = (im2SIFT_pts(1:2,matches(2,:)))';
 
matched_im1 = [matched_im1;inlier_im1.Location];  % total set of matches from sift and surf
matched_im2 = [matched_im2;inlier_im2.Location];

figure,showMatchedFeatures(im1,im2,inlier_im1,inlier_im2,'montage')
% showCorresponcesVertical(im1,im2,matched_im1,matched_im2,false)

tform = computeTransformation(inlier_im2.Location,inlier_im1.Location);
% tform = computeTransformation(matched_im2,matched_im1);

% warping using affine transformation
outputView = imref2d(size(im1));
newimage  = imwarp(I2,tform,'OutputView',outputView,'FillValues',255);
% label13.rectImage = newimage;


% % transforming the pts for the second image as well using tform
% T=tform.T';
% sz = size(label13.rectPoints{1});
% x=label13.rectPoints{1}(:);
% y=label13.rectPoints{2}(:);
% z=ones(length(x),1);
% oldcords = [x y z];
% newcords = round(T*oldcords');
% x=reshape(newcords(1,:),sz);
% y=reshape(newcords(2,:),sz);
% label13.rectPoints{1} = x;
% label13.rectPoints{2} = y;

end