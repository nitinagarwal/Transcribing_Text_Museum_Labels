function label = vertical_line_9(label)
% ver 8: trying nicola and marks suggestion.
% ver 9: updated as the textlines and WSlines are refined 

%---------------------
%vertical line...
%---------------------
nBin=180;

img=label.img_Correct;

[gx,gy]=imgradientxy(img(:,:,1));

% orientation gives the direction of the gradient, not the direction of text stroke. 
[mag, binning,ori] = im_bin(gx, gy, nBin); 

% both are perpendicular.

%removing the edges from magnitude. 
for k=1:4
    edg = [label.edges{k}(:,1) label.edges{k}(:,2)];
        for i=1:length(edg)
            mag(edg(i,2)-5:edg(i,2)+5,edg(i,1)-5:edg(i,1)+5)=0;
        end
end

if(size(mag)~=size(binning))
    error('Problem removing edges');
end

% dividing image into patches.

% lines=round(length(label.trajAll)/num_Y);
% ind=1:lines:length(label.trajAll);
% for i=1:length(ind)
%  pts(i,:,:)=interparc(num_X,label.trajAll{ind(i)}(:,1),label.trajAll{ind(i)}(:,2));
% end

num_X=15;
label.gridpts_perLine = num_X;

%chosing pts on the lines.
% z=1;
for i=1:length(label.ref_traj_textLines)
   
   vec = round(linspace(1,length(label.ref_traj_textLines{i}),num_X)); % sampling 
   pts(i,:,:) = label.ref_traj_textLines{i}(vec,:);
%    pts_inbound(z,:) = label.hori_Inbound{i}(vec,:);
   pts_blank(i,:) = label.pts_blank{i}(vec);
%    pts(z,:,:)=interparc(num_X,label.hori{i}(:,1),label.hori{i}(:,2));
%    z=z+1;
end

num_Y=size(pts,1);  % number of rows/textLines

%visualize
% figure,imshow(img,[]);hold on
% plot(pts(:,:,1)',pts(:,:,2)','g*');
% for i=1:num_Y
%     for j=1:num_X
%         if(pts_blank(i,j)==0)
%             plot(pts(i,j,1)',pts(i,j,2)','r*');
%          end
%     end
% end

step_X = round((pts(1,end,1)-pts(1,1,1))/num_X);
step_Y = round((pts(end,1,2)-pts(1,1,2))/size(pts,1)); %+10; % adding an additional constast to select the region more vertically

% step_X=step_X/2;
% step_Y=step_Y/2;

% histogram={};
% get the histogram around each point.
for i=1:num_Y
    for j=1:num_X
        
        patch_mag = mag(pts(i,j,2)-step_Y:pts(i,j,2)+step_Y,pts(i,j,1)-step_X:pts(i,j,1)+step_X);
        patch_bin = binning(pts(i,j,2)-step_Y:pts(i,j,2)+step_Y,pts(i,j,1)-step_X:pts(i,j,1)+step_X);
        patch_ori = ori(pts(i,j,2)-step_Y:pts(i,j,2)+step_Y,pts(i,j,1)-step_X:pts(i,j,1)+step_X);
        
         mag_histogram(:,i,j) = get_histogram2(patch_mag, patch_bin, nBin);
        [X,Y]=meshgrid(pts(i,j,1)-step_X:pts(i,j,1)+step_X,pts(i,j,2)-step_Y:pts(i,j,2)+step_Y);
    end
end


% computing the orientation of both the short edges.
% assumption no warping on the short edges hence fit linear line
[left_pts,left_err,left_coeff] = linearFit(label.edges{1});

left_angle = atand((left_coeff(1)));         % abs value cause slope can be -ve or +ve
% if(left_angle<0)                    % u need to flip it
%     left_angle = left_angle+180;
% end
left_angle = left_angle + 90;    % cause the direction is the gradient direction


[right_pts,right_err,right_coeff] = linearFit(label.edges{3});

right_angle = atand((right_coeff(1)));
% if(right_angle<0)                    % u need to flip it
%     right_angle = right_angle+180;
% end
right_angle = right_angle + 90;

% % if it is blank point, its direction is the mean of the direction of the
% % both the vertical edges
% 
% p=[label.edges{1}(1,:)' label.edges{1}(end,:)'];
% p=p(:,2)-p(:,1);
% left_angle=atand(p(2)/p(1));
% 
% p=[label.edges{3}(1,:)' label.edges{3}(end,:)'];
% p=p(:,2)-p(:,1);
% right_angle=atand(p(2)/p(1));  % taking abs as the slope can be -ve indicating going down

% mean vertical orientation
% mean_angle=mean([left_angle right_angle]);

% if(mean_angle < 0)
%     mean_angle = mean_angle + 180;
% end

% mean_angle = mean_angle + 90;  % cause the direction is the gradient direction

edges = linspace(0, 180, nBin + 1);
edgesBoundary = (edges(1:end-1) + edges(2:end)) / 2;

id=find(edgesBoundary<=left_angle);
leftDir=id(end);            
id=find(edgesBoundary<=right_angle);
rightDir=id(end);

for i = 1:num_Y  % # of textLines
     
        % nicola idea
       [dirs(i,:),weights(i,:)] = interpolatingEdges(mag_histogram(:,i,:),leftDir,rightDir,15,nBin,...
           pts_blank(i,:),pts(i,:,:));
        
        %show top 10 for each position
%         dir = topTENorientation(mag_histogram(:,i,:),leftDir,rightDir,15,nBin,pts_blank(i,:),...
%             pts(i,:,:));
% % %         
% % %         % marks idea (%take avg of top 10 values)
%         dir = dir-90;
%         dir(dir<0) = dir(dir<0)+180;
%         dir = mean(dir,1);
%         dirs = dir+90; 
%         
% %         %visualize top rank orientation
%         visual_vertical_cluster(pts(i,:,:),dirs);
end;

% there might be angle whose values are 180 apart (even though they point in same direction)
% . To solve those find >90 and subtract 180 from it. 
id = find(dirs >90);
dirs(id)=dirs(id)-180;


% making orientation same along small width
% hence taking a weighted average 
weights = 1./weights;   
weights(weights==Inf) = 100; % high weights
dirs = sum(weights.*dirs,1)./sum(weights,1);
dirs = repmat(dirs,num_Y,1);



% for i=1:num_Y
%     for j=1:num_X
%         if(pts_blank(i,j)==0)
%             dirs(i,j)=defaulDir;
%         end
%     end
% end

%visualize
% figure,imshow(img,[]);hold on
% for i=1:num_Y
%    visual_vertical(pts(i,:,:),dirs(i,:),false);
% %     pause
% end

% % averageing the vertical orientation. a single value per dir column.
% average_value = round(mean(dirs,1));
% dirs = repmat(average_value,num_Y,1);



for i=1:num_Y
    label.ver{i} = squeeze(pts(i,:,:));
    label.verOri{i}=dirs(i,:);
end

end
    










