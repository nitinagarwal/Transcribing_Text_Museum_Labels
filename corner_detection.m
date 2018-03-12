function label = corner_detection(label)
%Input - C is the input image.
%Output - left and right corner points

C = label.img;
debug=label.debug;

% correcting for rotation using resampled convex hull and computing the
% convex hull of the rotated image
%  [Icorrected,theta]=rotationCorrection1(C,debug);
Icorrected = C;

label.mask = maskComputation(Icorrected);
label.img_Correct=Icorrected;

rotated_convexhull = rotated_convexHull(label.mask ,debug);

convexhullImage = zeros(size(label.mask));
for i=1:size(rotated_convexhull,2)
    convexhullImage(round(rotated_convexhull(2,i)),round(rotated_convexhull(1,i))) = 1;
end

se = strel('disk',3);
convexhullImage = imdilate(convexhullImage,se);
convexhullImage = imfill(convexhullImage);
se = strel('disk',5);
convexhullImage = imerode(convexhullImage,se);
[B,~,~,~] = bwboundaries(convexhullImage,'noholes');

%these are the boundary points
boundary=B{1};
convexhullImage = zeros(size(label.mask));
id=sub2ind(size(convexhullImage),boundary(:,1),boundary(:,2));
convexhullImage(id)=1;
label.convexhull = convexhullImage;

gray=rgb2gray(Icorrected);
edgeImage=edge(gray,'canny'); % let ostu's method determine the threshold
% figure,imshow(edgeImage)

% removing small debri
cleanImage=largestConnectedComponent(edgeImage,50,debug);

% modified Hough Transform
lines = modified_Hough_Transform(cleanImage,debug);

if(debug)
    figure,imshow(gray);hold on
    for k = 1:length(lines)
       xy = [lines(k).point1; lines(k).point2];
       plot(xy(:,1),xy(:,2),'LineWidth',1,'Color','green');

       % Plot beginnings and ends of lines
       plot(xy(1,1),xy(1,2),'mx','LineWidth',2);
       plot(xy(2,1),xy(2,2),'mx','LineWidth',2);
    %    pause;
    end
end

% Computing 2 lines representing 2 edges of the label.(vertical edges)
% vertical lines less likey to bend. hence hough transform would most
% likely work.

[line1,line2]=vertical_lines(lines,Icorrected,debug);

% plot final lines
figure,imshow(Icorrected);hold on

l1 = line_equation(line1(1,:),line1(2,:));
pts1 = [-(l1(2)+l1(3))/l1(1) 1 ; -(l1(3)+(l1(2)*size(gray,1)))/l1(1) size(gray,1)];
plot(pts1(:, 1), pts1(:, 2), '-g','LineWidth',2);

l2 = line_equation(line2(1,:),line2(2,:));
pts2 = [-(l2(2)+l2(3))/l2(1) 1 ; -(l2(3)+(l2(2)*size(gray,1)))/l2(1) size(gray,1)];
plot(pts2(:, 1), pts2(:, 2), '-g','LineWidth',2);

% plotting the small lines 
plot(line1(:,1),line1(:,2),'-r','LineWidth',2);
plot(line2(:,1),line2(:,2),'-r','LineWidth',2);


%%%%%%%%%%%% horizontal lines more likely to bend hence hough transform
%%%%%%%%%%%% less likely to work. Therefore find corners via intersection of
%%%%%%%%%%%% convex hull and the two veretical lines.

left_points = horizontal_lines(line1,line2,rotated_convexhull,'left');
right_points = horizontal_lines(line1,line2,rotated_convexhull,'right');

label.cornerpts{1}=left_points;
label.cornerpts{2}=right_points;

plot(left_points(:,1),left_points(:,2),'*b','MarkerSize',10)
plot(right_points(:,1),right_points(:,2),'*b','MarkerSize',10)

l3 = line_equation(left_points(1,:),right_points(1,:));
pts3 = [-(l3(2)+l3(3))/l3(1) 1 ; -(l3(3)+(l3(2)*size(gray,1)))/l3(1) size(gray,1)];
plot(pts3(:, 1), pts3(:, 2), '-g','LineWidth',2);

l4 = line_equation(left_points(2,:),right_points(2,:));
pts4 = [-(l4(2)+l4(3))/l4(1) 1 ; -(l4(3)+(l4(2)*size(gray,1)))/l4(1) size(gray,1)];
plot(pts4(:, 1), pts4(:, 2), '-g','LineWidth',2);

disp('Left Corner Points') ;
disp(left_points);
disp('Right Corner Points');
disp(right_points); 

end


