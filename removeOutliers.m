function PTS = removeOutliers(in_pts)
% removing outliers based on tangent to the left and right neighbour

ycords = in_pts(:,2);

next = circshift(ycords,-1);
prev = circshift(ycords,1);

error = [abs(ycords-next) abs(ycords-prev)];

threshold = 3.0;

id = find(error(:,1) <= threshold & error(:,2) <= threshold);

PTS = in_pts(id,:);

end