function [A,xy]=img2graph_modified(B,edgeImage)
% A is the adjacency matrix of image. This is modified according to the
% edgeImage. 
% xy is the coordinates of all the vertices.
% B is a matrix which stores the position of the vertices.
% edgeImage is the image whose connectivity needs to be applied to A.


[M,N] = size(B);                      %# Get the matrix size

A = sparse(M*N,M*N);   % creating a sparse matrix with zeros

n_zeros=B.*edgeImage;

values=nonzeros(n_zeros);

% padding to avoid boundary cases.
n_zeros=padarray(n_zeros,[1 1],0);

% creating the adjaceny matrix
for i=1:length(values)
    
   [col,row]=ind2sub([size(B,2) size(B,1)],values(i)); % this needs to be reversed(check B)
   
   temp=n_zeros((row+1)-1:(row+1)+1,(col+1)-1:(col+1)+1); % because of padding.
   temp=nonzeros(temp);
   [id,~]=find(temp==values(i)); % removing self;
   temp(id)=[];
   
   A(values(i),temp)=1;
   A(temp,values(i))=1;
end

%creating the cost matrix
% cost = ones(M*N,M*N);


% creating the cooridates array.
[X,Y]=meshgrid(1:N,1:M);
x=reshape(X',[],1);
y=reshape(Y',[],1);
xy=[x y];


end

