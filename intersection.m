function   intersect_point = intersection(points1,points2)
% Intersection of two small line segments.
% Input is 2 point sets having two points each. total 4 points
% points row-wise
% points1 = from the straight vertical line
% points2 = from the convex hull

l1 = line_equation(points1(1,:),points1(2,:));
l2 = line_equation(points2(1,:),points2(2,:));

%plotting the line
% pts = [-(l2(2)+l2(3))/l2(1) 1 ; -(l2(3)+(l2(2)*size(gray,1)))/l2(1) size(gray,1)];
% plot(pts(:, 1), pts(:, 2), '-m','LineWidth',1);

inter=cross(l1,l2);

point_x=inter(1)/inter(3);
point_y=inter(2)/inter(3);

%checking whether the point lies in between the two line segments (both vectors are 180 apart).

vec1 = [point_x-points2(1,1) point_y-points2(1,2)];
vec2 = [point_x-points2(2,1) point_y-points2(2,2)];

result = dot(vec1,vec2)/(norm(vec1)*norm(vec2));   % cos(theta)

if(round(result) == -1 || norm(vec1)==0 || norm(vec2)==0 )   % round to remove any numerical error.
    intersect_point = [point_x point_y];
else
    intersect_point = [nan nan];
end

% % if (point_x >= points1(1,1) && point_x <= points1(2,1) && point_y >= points1(1,2) && point_y >= points1(2,2))
%     % lies on first segment
%     if (point_x >= points2(1,1) && point_x <= points2(2,1) && point_y >= points2(1,2) && point_y >= points2(2,2) )
%         
%         intersect_point = [point_x point_y];
%     else
%         intersect_point = [0 0];
%     end
%     
% % else
% %     intersect_point = [0 0];
% % 
% % end

end