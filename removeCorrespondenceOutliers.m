function [inlier1,inlier2] = removeCorrespondenceOutliers(pts1,pts2,im1,im2)
% input: first pt- pts1
%      : second pt-  pts2
%output: [inlier1 inlier2]

radius = 10;  %was 100
numCorr = size(pts1,1);
inlier1=[];
inlier2=[];

h1 = size(im1,1);
h2 = size(im2,1);
h = min(h1,h2);

[X1,~]=meshgrid(1:size(im1,2),1:size(im1,1));
[X2,~]=meshgrid(1:size(im2,2),1:size(im2,1));

for i=1:numCorr
  
    %template img in 1
    filter1 = im1( 1:h, X1(pts1(i,2),pts1(i,1))-radius : X1(pts1(i,2),pts1(i,1))+radius );
    filter2 = im2( 1:h, X2(pts2(i,2),pts2(i,1))-radius : X2(pts2(i,2),pts2(i,1))+radius );
    score1 = computeDiff(filter1, filter2); % score with correspondence
    C=normxcorr2(filter1,filter2);
    score3 = max(C(:));
    if(score3 > 0.7)
        disp('score > 0.7');
    end
    
    filter2 = im2( 1:h, X1(pts1(i,2),pts1(i,1))-radius : X1(pts1(i,2),pts1(i,1))+radius );
    score2 = computeDiff(filter1, filter2); % score before correspondence
    C=normxcorr2(filter1,filter2);
    score4 = max(C(:));
    
%     if ( score1-score2 >0.01 || score1 > 0.6) 
    if ( score3 > 0.7)     % correspondence is making it score worse
%        score2/score1 >=1.2 || 
        inlier1=[ inlier1; pts1(i,:)];
        inlier2=[ inlier2; pts2(i,:)];
        
    end

end
% less score is better - more similar the two regions are
end

function diff = computeDiff(im1, im2)
%input are grayscale images

assert(size(im1,1)==size(im2,1));
assert(size(im1,2)==size(im2,2));

diff = (im1-im2).^2;
diff = sqrt(sum(diff(:)));

end


function NCC_score = normalizedCross(im1, im2)
% im1 - fixed image
% im2 - moving image

% check both are same dimns and gray images
assert(ndims(im1)==ndims(im2));
assert(ndims(im1)~=3);

a=im2double(im1);
template = im2double(im2);

% bluring template
havg=fspecial('average',[2,2]);
template=conv2(template,havg,'same');

% normalized correlation
C=normxcorr2(template,a); 
% figure,imagesc(C)
% colormap(jet)

% % thresholding.
% C(C<=0.54)=0;
% 
% %nonmaximal suppression
% for i=2:size(C,1)-1
%     for j=2:size(C,2)-1
%       
%         if(C(i,j)<=C(i-1,j-1) | C(i,j)<=C(i,j-1) | C(i,j)<=C(i+1,j-1)...
%                 |C(i,j)<=C(i-1,j) |C(i,j)<=C(i+1,j) |C(i,j)<=C(i-1,j+1)...
%                 |C(i,j)<=C(i,j+1) |C(i,j)<=C(i+1,j+1))
%             C(i,j)=0;
%         end
%     end
% end

% detecting the location of max value.

NCC_score = max(C(:));

% [r,c]=find(C==max(C(:)));

% r_offset = r - size(template,1);
% c_offset = c - size(template,2);


% r=r-size(template,1)/2-0.5;
% c=c-size(template,2)/2-0.5;


% figure,imshow(a);hold on
% plot(c,r,'r*')
end