function label = edge_detectionNew(label)
% Computing the edge pixels from the four corner points
% No need to use the dijksta algorithm.

flag=label.debug;
Icorrected = label.img_Correct;

% rounding the corner pts.
left_points1=round(label.cornerpts{1});
right_points1=round(label.cornerpts{2});

% edgeImage = label.mask;     
edgeImage = label.convexhull;   % just for wacv2018

% insuring that cornerpts lie on the boundary
[left_points1,right_points1] = modify_edgeImage(edgeImage,left_points1,right_points1);

% storing modified cornerpts
label.cornerpts{1}=left_points1;
label.cornerpts{2}=right_points1;

if(flag)
    figure,imshow(edgeImage);hold on
    plot(left_points1(:,1),left_points1(:,2),'*b','MarkerSize',10)
    plot(right_points1(:,1),right_points1(:,2),'*b','MarkerSize',10)
end

% order is: left_top;left_bottom;right_top;right_bottom
points = [left_points1(1,:);left_points1(2,:);right_points1(1,:);right_points1(2,:)];

edgels{1} = wavy_edgeDetection(edgeImage,points(1,:),points(2,:)); % left
edgels{2} = wavy_edgeDetection(edgeImage,points(3,:),points(1,:)); % top
edgels{3} = wavy_edgeDetection(edgeImage,points(4,:),points(3,:)); % right
edgels{4} = wavy_edgeDetection(edgeImage,points(4,:),points(2,:)); % bottom

% checking whether path was found or not.
if(isnan(edgels{1}(1,1)) || isnan(edgels{2}(1,1)) || isnan(edgels{3}(1,1)) || isnan(edgels{4}(1,1)) )
    error ('Edges are not computed properly as they are disconnected');
end

figure,imshow(Icorrected);hold on
color=['r';'g';'b';'m'];
for k=1:4
plot(edgels{k}(:,1),edgels{k}(:,2),[color(k) '*'],'MarkerSize',1);
end

label.edges=edgels;

% flip the left edge as its computed in reverse
label.edges{1}=flipud(label.edges{1});

end









