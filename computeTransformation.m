function tform = computeTransformation(pts1,pts2)
% Input:
% are the correspondences
% Output:
% transformation affine2d matrix


A = [pts1 ones(length(pts1),1)];
bx = pts2(:,1);
by = pts2(:,2);

[U,S,V]=svd(A,0);
q = V*(S\(U'*bx));
r = V*(S\(U'*by));

transformation=[q';r'];

tform=affine2d(transformation');

end