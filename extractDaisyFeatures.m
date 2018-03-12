function features = extractDaisyFeatures(dzy,points)
%output is a feature vector with valid pts .

% locations = points.Location();
locations = round(points);
locations = unique(locations,'rows');

x=locations(:,1);
y=locations(:,2);

% function out=display_descriptor(dzy, y, x)
z=1;
for i=1:length(x)
    
if y(i)<0 || x(i)<0 || y(i)>dzy.h-1 || x(i)>dzy.w-1
    continue
else

    validPoints(z,:) = [x(i) y(i)];   % valid pts
    out = reshape( dzy.descs( y(i)*dzy.w+x(i)+1, :), dzy.HQ, dzy.HN )';
    featureVector(z,:)=reshape(out,[],1);
    z=z+1;
end

features.vector = featureVector;
features.location = validPoints;


end