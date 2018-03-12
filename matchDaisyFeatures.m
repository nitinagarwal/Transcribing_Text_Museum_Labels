function matchedpts = matchDaisyFeatures(featuresDaisy_im1,featuresDaisy_im2,dzy1,dzy2,im1,im2)

radius = 20; % was 20
threshold = 1.0;
matchpts_im1=[];
matchpts_im2=[];
score = [];

fvector1 = featuresDaisy_im1.vector;
fvector2 = featuresDaisy_im2.vector;
l1 = featuresDaisy_im1.location;
l2 = featuresDaisy_im2.location;

for i=1:size(fvector1,1)
    
    if(l1(i,1)<radius+1 || l1(i,2)<radius+1)
        continue;
    end
    
    [X,Y]=meshgrid(l1(i,1)-radius:l1(i,1)+radius,l1(i,2)-radius:l1(i,2)+radius); % computing via neighbourhood
    p=[X(:) Y(:)];
    features = extractDaisyFeatures(dzy2,p);
    result1 = bsxfun(@minus,features.vector,fvector1(i,:));  %eucledian distance between feature vector
    result1 = sqrt(sum(result1.^2,2));
    id1 = find(result1==min(result1));
    
    if(numel(id1)>1)
        disp('id1');
%         break;
        [id1,~] = findBetterNCC(l1(i,:),features.location(id1,:),im1,im2);
%         [id1,~] = findClosestPoint(l1(i,:),features.location(id1,:));       % closest pt if tie
    end
    
    result2 = bsxfun(@minus,fvector2,fvector1(i,:));  %eucledian distance between feature vector
    result2 = sqrt(sum(result2.^2,2));
    id2 = find(result2==min(result2));
    
    if(numel(id2)>1)
        disp('id2');
%         break;
        [id2,~] = findBetterNCC(l1(i,:),features.location(id1,:),im1,im2);
%         [id2,~] = findClosestPoint(l1(i,:),featuresDaisy_im2.location(id2,:)); % closest pt if tie
    end
    
    if(result1(id1) < result2(id2))
       
            matchpts_im1 = [matchpts_im1; l1(i,:)];
            matchpts_im2 = [matchpts_im2;features.location(id1,:)];
            score = [score;result1(id1)];
        
    elseif(result1(id1) > result2(id2))
    
        matchpts_im1 = [matchpts_im1; l1(i,:)];
        matchpts_im2 = [matchpts_im2;featuresDaisy_im2.location(id2,:)];
        score = [score;result2(id2)];
        
    elseif(result1(id1) == result2(id2))
            
           disp('score is equal');
%            break;
           [value,~] = findBetterNCC(l1(i,:),[features.location(id1,:); featuresDaisy_im2.location(id2,:)],im1,im2);
%          [value,~] = findClosestPoint(l1(i,:),[features.location(id1,:); featuresDaisy_im2.location(id2,:)]);
            
           if(value==1)
               matchpts_im1 = [matchpts_im1; l1(i,:)];
               matchpts_im2 = [matchpts_im2;features.location(id1,:)];
           else
               matchpts_im1 = [matchpts_im1; l1(i,:)];
               matchpts_im2 = [matchpts_im2;featuresDaisy_im2.location(id2,:)];
           end
           score = [score;result1(id1)];
    end
%     pause
%     figure,imshow(im1);hold on
%     plot(l1(i,1),l1(i,2),'b*','MarkerSize',3);
%     figure,imshow(im2);hold on
%     plot(matchpts_im2(end,1),matchpts_im2(end,2),'r*','MarkerSize',3);
    
end

if(size(matchpts_im1,1) <5)
    matchedpts.score = nan;
    matchedpts.image1=nan;
    matchedpts.image2=nan;
    return;
end


% outlier removal (ransac)
% assuming there is an affine transformation between two images using just
% these small correspondences.
% [~, inlier_im1,inlier_im2] = estimateGeometricTransform(matchpts_im1,matchpts_im2,'affine','MaxDistance',20);

% sprintf('size of inliers is %i',size(inlier_im1,1))
% [inlier_im1,inlier_im2] = removeCorrespondenceOutliers(inlier_im1,inlier_im2,im1,im2);
[inlier_im1,inlier_im2] = removeCorrespondenceOutliers(matchpts_im1,matchpts_im2,im1,im2);
sprintf('size of modified inliers is %i',size(inlier_im1,1))

if(isempty(inlier_im1))
    matchedpts.score = nan;
    matchedpts.image1=nan;
    matchedpts.image2=nan;
    return;
end

[~,~,ib]=intersect(inlier_im1,matchpts_im1,'rows');
newscore = score(ib);

% sprintf('size of inliers is %i',length(newscore));
% try distribute the pts evenly or randomly
% numpts = round(length(newscore)/10);
numpts=1;
matchedpts.score = newscore(1:numpts:end);
matchedpts.image1=inlier_im1(1:numpts:end,:);
matchedpts.image2=inlier_im2(1:numpts:end,:);


% attempt 1
%taking 10 pts with the lowest scores
% [score,id]=sort(newscore,'descend');

% matchedpts.score = score(1:10);
% matchedpts.image1=inlier_im1(id(1:10),:);
% matchedpts.image2=inlier_im2(id(1:10),:);

% attempt 2
% %based on distance
% dis = sqrt(sum((matchpts_im1 - matchpts_im2).^2,2));
% 
% id = find(dis<100);
% 
% a = matchpts_im1(id,:);
% b = matchpts_im2(id,:);




end