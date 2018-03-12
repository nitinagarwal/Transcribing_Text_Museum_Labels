function  pts= estimate_perspective(pt1,pt2,query)
% Computing the projective transformation matrix for these corresponding pts (DLT algorithm)
% pt1 -> pt2

% id=1;
% for i=1:length(pt1)
%     
%   A(id,:) = [-pt1(i,1) -pt1(i,2) -1 0 0 0 pt1(i,1)*pt2(i,1) pt1(i,2)*pt2(i,1) pt2(i,1)];
%   A(id+1,:) = [0 0 0 -pt1(i,1) -pt1(i,2) -1  pt1(i,1)*pt2(i,2) pt1(i,2)*pt2(i,2) pt2(i,2)];
%   id=id+2;
% end
% 
% % computing the transformation matrix T
% [~,~,V]=svd(A,0);
% H=V(:,9);
% T= reshape(H,[3 3])';


T = DirectLinearTransformation(pt1',pt2');

query = padarray(query',[1 0],1,'post');

temp_pts=T*query;

% homogenous to cartesian
pts(1,:) = temp_pts(1,:)./temp_pts(3,:);
pts(2,:) = temp_pts(2,:)./temp_pts(3,:);

pts=pts';




end





