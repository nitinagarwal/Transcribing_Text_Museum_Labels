 function points = horizontal_lines(line1,line2,rotated_convexhull,flag)
% line1 is the left vertical line points. line2 is the right vertical line
% points.
% rotated_convexhull is the convex hull of the rotated image. 
% flag is whether you want left corner points or right corner points.

new1=rotated_convexhull;
new1(:,end+1)=new1(:,1);  % circular array


if (strcmp(flag,'left'))

    for i=1:size(new1,2)-1
        intersect_point(i,:) = intersection(line1,[new1(:,i)';new1(:,i+1)']);        
    end 
    
elseif(strcmp(flag,'right'))
    
    for i=1:size(new1,2)-1
        intersect_point(i,:) = intersection(line2,[new1(:,i)';new1(:,i+1)']);        
    end
    
else 
    error('Which corner points you want');
end


% getting the top most and bottom most intersect_points and returning the
% points in order. (works for both left and right points)

id = find(~isnan(intersect_point(:,1)));
if(numel(id)>2 )
   
 low_id = find(intersect_point(id,2) == min(intersect_point(id,2)) );
 high_id = find(intersect_point(id,2) == max(intersect_point(id,2)) );
    
   points=[intersect_point(id(low_id),:);intersect_point(id(high_id),:)];
      
else 
    points=[intersect_point(id(1),:);intersect_point(id(2),:)];
end

% flipping it, if found in reverse
if(points(2,2)<points(1,2))
   points=flipud(points);
end


    

 end

 
 
 
 
 