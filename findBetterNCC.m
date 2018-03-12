function [id,dis] = findBetterNCC(pt,otherpts,im1,im2)
% input: first pt- target pt
%      : second pt- query pts
%output: which query pt and its NCC score

assert(size(otherpts,1)>1,'comparison pts less than 2');
radius = 100; % was 100

h1 = size(im1,1);
h2 = size(im2,1);
h = min(h1,h2);

[X1,Y1]=meshgrid(1:size(im1,2),1:size(im1,1));
[X2,Y2]=meshgrid(1:size(im2,2),1:size(im2,1));

%template img in 1
filter1 = im1( 1:h, X1(pt(2),pt(1))-radius(1):X1(pt(2),pt(1))+radius(1) );
% filter1 = im1( Y1(pt(2),pt(1))-radius(2):Y1(pt(2),pt(1))+radius(2), X1(pt(2),pt(1))-radius(1):X1(pt(2),pt(1))+radius(1) );

for i=1:length(otherpts)
   
    filter2 = im2(1:h, X2(otherpts(i,2),otherpts(i,1))-radius(1):X2(otherpts(i,2),otherpts(i,1))+radius(1) );
    score(i) = computeDiff(filter1, filter2);
    
end
% less score is better - more similar the two regions are

[dis,id]=sort(score,'ascend');  %shortest distance
id = id(1);
dis = dis(1);



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