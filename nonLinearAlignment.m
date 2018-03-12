function  [label1,label2] = nonLinearAlignment(label1,label2)

% reference image is the one with more number of textLines
if(length(label1.ver)>=length(label2.ver))
    R = label1;
    M = label2;
    flag = true;
else
    R = label2;
    M = label1;
    flag = false;
end

% assert(size(R,1)==size(M,1) && size(R,2)==size(M,2));

textLines1 = median(R.rectPoints{2},2);
textLines2 = median(M.rectPoints{2},2);
P1=[]; % correspondences for the whole image
P2=[];

tic
for j=1:2:length(textLines2)   

im1 = R.affine(textLines1(j)-50:textLines1(j+1)+50,:,:);    % adding some buffer 
im2 = M.affine(textLines2(j)-50:textLines2(j+1)+50,:,:);    
im1=rgb2gray(im1);
im2=rgb2gray(im2);


%--------------------------feature detection & extraction-----------
p1=[];
p2=[];
im1_pts = detectMSERFeatures(im1);
im2_pts = detectMSERFeatures(im2);

p1=[p1;im1_pts.Location];
p2=[p2;im2_pts.Location];

im1_pts = detectHarrisFeatures(im1,'MinQuality',.01);
im2_pts = detectHarrisFeatures(im2,'MinQuality',.01);

p1=[p1;im1_pts.Location];
p2=[p2;im2_pts.Location];

p1 = round(p1);
p1 = unique(p1,'rows');         % removing any duplicates if there are due to rounding off
p2 = round(p2);
p2 = unique(p2,'rows');


dzy1 = compute_daisy(im1);
dzy2 = compute_daisy(im2);

%---------------------------------------dividing the image into cells
% both the images will have width same as they are affinely aligned. 

%cell orientations: 
% 1  2  3  4  5  6  7  8  9  10
% 11 12 13 14 15 16 17 18 19 20

cell1 = divideintoCells(im1);   % always 20 cells 
cell2 = divideintoCells(im1);
pts1 = [];             %pts in image 1
pts2 = [];             %pts in image 2

for i=1:length(cell1)
    %for every cell check first if its has any text and any harris corners. 
    % then compute features for only those corner location and do matching
    % with the other 3 neighbouring cells in the other image
    if(cell1(i).isText)
       
        validpts = intersect(cell1(i).Location,p1,'rows');
        if(size(validpts,1)==0)
            continue;
        end
        
        featuresDaisy_im1 = extractDaisyFeatures(dzy1,validpts);
        
        if(i<=10)
           % neighbours are down
           featuresDaisy_im2 = extractDaisyFeaturesNeighbours(cell2,i,p2,dzy2,true);
        else
          % neighbours are up
           featuresDaisy_im2 = extractDaisyFeaturesNeighbours(cell2,i,p2,dzy2,false);
        end
        
        if(size(featuresDaisy_im1.location,1)==0 || size(featuresDaisy_im2.location,1)==0 )
            continue;                           % no pts to match
        end
        
       indexPairs = matchDaisyFeatures(featuresDaisy_im1,featuresDaisy_im2,dzy1,dzy2,im1,im2); 
       if(~isnan(indexPairs.score))
           pts1 = [pts1;indexPairs.image1];
           pts2 = [pts2;indexPairs.image2];
       end
    end
end

if(~isempty(pts1))      % only we found some points
    P1=[P1; [pts1(:,1) pts1(:,2)+textLines1(j)-50] ];
    P2=[P2; [pts2(:,1) pts2(:,2)+textLines2(j)-50] ];   
end

end

% remove duplicates from P1 & P2 separately if any (regions overlapped and
% hence there will some be duplicates)
[P1,id] = unique(P1,'rows');
P2=P2(id,:);

[P2,id] = unique(P2,'rows');
P1=P1(id,:);

% outlier detection. (Assump: there are few outliers)
delta_x = P1(:,1)-P2(:,1);
delta_y = P1(:,2)-P2(:,2);

fit_x = medfilt1(delta_x,50,'omitnan','truncate');
fit_x = abs(fit_x-delta_x);
id1 = find(fit_x>25);                   %change allowed is 25          

fit_y = medfilt1(delta_y,50,'omitnan','truncate'); 
fit_y = abs(fit_y-delta_y);
id2 = find(fit_y>10);                   %change allowed is 10  (less for y)

% up = mean(delta_x)+std(delta_x);
% dwn= mean(delta_x)-std(delta_x);
% 
% id1 = find(delta_x>up | delta_x<dwn);
% up = mean(delta_y)+10;                  % random variable 10
% dwn= mean(delta_y)-10;
% id2 = find(delta_y>up | delta_y<dwn);
id=[id1;id2];                           % combining the outliers
id=unique(id);                          % there might be some duplicates

P1(id,:)=[];
P2(id,:)=[];

% take correspondences uniformly across image

showCorresponcesVertical(R.affine,M.affine,P1,P2,false)

% final non-linear warping
new_image = LaplaceWarping(P2,P1,M.affine,R.affine);
R.aligned = R.affine;
M.aligned = new_image;
figure,imshowpair(R.aligned ,M.aligned,'falsecolor');

fprintf('Time for Alignment is %f secs \n',toc) 

% reference image is the one with more number of textLines
if(flag)
    label1 = R;
    label2 = M;
else
    label2 = R;
    label1 = M;
end


end