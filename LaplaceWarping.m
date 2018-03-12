function new_image = LaplaceWarping(pt1,pt2,C1,C2)
% Input:
% pts1 are the initial points.
% pts2 are the final points where we would like them to be moved
% C1 is the input_image which is to be warped.
% C2 is the input_image which is the reference

% Output:
% outputimage: warped output image.
% For interpolation we are computing a delauny traingulation and then
% interpolating(refer to scatteredinterpolation documentation in matlab)

pt1 = round(pt1);
pt2 = round(pt2);

A=[pt1 ones(size(pt1,1),1)];
bx=pt2(:,1);
by=pt2(:,2);
I=C1;

% creating three images. boundary, xdis, ydis
boundary=zeros(size(I,1),size(I,2));                            

for index=1:length(A)             % input image having pixel value 1 where we have correspondece
    boundary(A(index,2),A(index,1))=1;
end

xdis=zeros(size(I,1),size(I,2)); 
ydis=zeros(size(I,1),size(I,2));

for index=1:length(bx)
    
    xdis(A(index,2),A(index,1))= bx(index) - A(index,1);        % x and y displacement images
    
    ydis(A(index,2),A(index,1))= by(index) - A(index,2);
end

[Lx, Ly]=solveLaplace(xdis,ydis,int32(boundary));               % Laplace Transform

x=zeros(size(I(:,:,1)));
y=zeros(size(I(:,:,1)));

for i=1:size(Lx,1)                           % rows of the image
    for j=1:size(Lx,2)                       % columns of the image
               
        y(i,j)= i+Ly(i,j) ;
        x(i,j)= j+Lx(i,j) ;

    end
end

%deformation for output image
[r,c,~]=size(C2);
[xq,yq]=meshgrid(1:c,1:r);

x=reshape(x,[],1);
y=reshape(y,[],1);

% tic
parfor i=1:3
% output_image(:,:,i) = griddata(x,y,double(I(:,:,i)),xq,yq,'linear');
F = scatteredInterpolant(x,y,reshape(double(I(:,:,i)),[],1),'linear','linear');
new_image(:,:,i)=F(xq,yq);
end
% fprintf('Time for laplace warping is %f secs \n',toc) 

new_image=uint8(new_image);

end

