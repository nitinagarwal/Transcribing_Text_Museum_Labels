function edgeImage = maskComputation(image)
%computes the BW mask of the label

gray=rgb2gray(image);
Inew=zeros(size(gray,1),size(gray,2));

id1 = find(gray<=240); % the label is already segmented hence no need to run canny
Inew(id1)=1;

Inew=largestConnectedComponent(Inew,10,false); % removing any debris


% edgeImage=edge(gray,'canny');
% Inew=largestConnectedComponent(edgeImage,100,flag);  % remove any debris if anything comes in background%%% CHECK

[B,~,~,~] = bwboundaries(Inew,'noholes');

%these are the boundary points
boundary=B{1};

edgeImage=zeros(size(gray,1),size(gray,2));
id=sub2ind(size(gray),boundary(:,1),boundary(:,2));
edgeImage(id)=1;

end