function [output_pts,bestModel] = ransac(input_pts)
% Input:
% input_pts ? a set of observed data points
% model ? a model that can be fitted to data points (cubic) -
% Ax^3+Bx^2+Cx+D.
% n ? the minimum number of data values required to fit the model
% k ? the maximum number of iterations allowed in the algorithm
% t ? a threshold value for determining when a data point fits a model
% d ? the number of close data values required to assert that a model fits well to data
% Output:
%  the number of inliers and the model parameters

if(size(input_pts,2)~=2)
    input_pts=input_pts';
end

num = size(input_pts,1);        % Total number of points

n = 4;
k = 5000;
t = 5;
d = round(0.90*num);            % searching for 90% of the input
    
bestModel = [];
bestInlierNum = 0;
besterr = exp(10);              % some very high value
output_pts = [];


for i=1:k

% Randomly select n points
     idx = randperm(num,n); sample = input_pts(idx,:);      
     model = polyfit(sample(:,1),sample(:,2),3); %fit cubic model
     other_idx = setdiff([1:num],idx);
    
%      potential_inlier = [];
% error of remaining pts to this model
        yValue =  polyval(model,input_pts(other_idx,1)); 
        dis = sqrt((input_pts(other_idx,2)-yValue).^2);
        potential_inlier = input_pts(other_idx(dis<t),:) ;

%     for j=1:length(other_idx)
%         yValue =  polyval(model,input_pts(other_idx(j),1)); 
%         dis = sqrt((input_pts(other_idx(j),2)-yValue).^2);
%         
%         if(dis<t)
%         potential_inlier = [potential_inlier; input_pts(other_idx(j),:)];   
%         end
%     end
           
%     figure,plot(input_pts(:,1),input_pts(:,2),'ro','MarkerSize',3);hold on
%     plot(sample(:,1),sample(:,2),'g*','MarkerSize',4);hold on

% Update the number of inliers and fitting model if better model is found  
    if(length(potential_inlier)>d)
        
        pts = [sample;potential_inlier];
        yValue =  polyval(model,pts(:,1)); % model parameters fitted to all points in maybeinliers and alsoinliers
        dis = sum(sqrt((pts(:,2)-yValue).^2));
%         disp('inlier length greater than d');
        
        if(dis < besterr)
            bestModel = model;
            besterr = dis;
%             bestInlierNum = length(potential_inlier);
            output_pts = pts;
            
%             disp('found a good model');  
%            plot(output_pts(:,1),polyval(bestModel,output_pts(:,1)),'b-','LineWidth',1)
        end
    end
    
%     pause
end

if(size(output_pts,1)==0)
    output_pts = input_pts;
    bestModel = [0 0];
    return 
end

%sorting the outputpts by x
[~,id] = sort(output_pts(:,1));
output_pts = output_pts(id,:);

% figure,plot(input_pts(:,1),input_pts(:,2),'ro','MarkerSize',3);hold on
% plot(output_pts(:,1),output_pts(:,2),'g*','MarkerSize',4);hold on
% plot(output_pts(:,1),polyval(bestModel,output_pts(:,1)),'b-','LineWidth',1)
% set(gca,'YDir','reverse')


%---------------------------------taken from wiki--------------------------
% iterations = 0
% bestfit = nul
% besterr = something really large
% while iterations < k {
%     maybeinliers = n randomly selected values from data
%     maybemodel = model parameters fitted to maybeinliers
%     alsoinliers = empty set
%     for every point in data not in maybeinliers {
%         if point fits maybemodel with an error smaller than t
%              add point to alsoinliers
%     }
%     if the number of elements in alsoinliers is > d {
%         % this implies that we may have found a good model
%         % now test how good it is
%         bettermodel = model parameters fitted to all points in maybeinliers and alsoinliers
%         thiserr = a measure of how well model fits these points
%         if thiserr < besterr {
%             bestfit = bettermodel
%             besterr = thiserr
%         }
%     }
%     increment iterations
% }
% return bestfit













end