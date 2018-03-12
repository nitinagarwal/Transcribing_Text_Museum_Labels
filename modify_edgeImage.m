function [left_points1,right_points1] = modify_edgeImage(edgeImage,left_points1,right_points1)
% with this implementation of Dijkstras algorithm this needs to be done.
% the corner points might not lie on edgeImage. Hence add those to the
% edgeImage. 

if( ~edgeImage(left_points1(1,2),left_points1(1,1)) )
    
    left_points1(1,:) = shiftPoint(edgeImage,left_points1(1,:));
end

if( ~edgeImage(left_points1(2,2),left_points1(2,1)) )
    
    left_points1(2,:) = shiftPoint(edgeImage,left_points1(2,:));
end

if( ~edgeImage(right_points1(1,2),right_points1(1,1)) )
    
    right_points1(1,:) = shiftPoint(edgeImage,right_points1(1,:));
end

if( ~edgeImage(right_points1(2,2),right_points1(2,1)) )
    
    right_points1(2,:) = shiftPoint(edgeImage,right_points1(2,:));
end

end

function Output_pt = shiftPoint(edgeImage,input_pt)

[r,c]=find(edgeImage==1);
p=[c r];

% computing distance of input_pt to all the points in edgeImage & chose min
dis = bsxfun(@plus,dot(p,p,2),dot(input_pt,input_pt,2))-2*p*input_pt';
id = find(dis==min(dis));

Output_pt=[p(id(1),1) p(id(1),2)];

end


