function convex_hull = rotated_convexHull(image,flag)
% given a image mask it computes the convex hull of the rotated image.


[rows1,cols1]=find(image==1);
k1=convhull(cols1,rows1);                           % computing the convex hull
new1=interparc(1000,cols1(k1),rows1(k1),'linear');  % resampling the convex hull at equal intervals (1000 points)

if(flag)
    figure,imshow(image);hold on;
    plot(new1(:,1),new1(:,2),'g*','MarkerSize',4); % resampled convex hull
end

convex_hull= new1';

end