function imageCell = divideintoCells(im)
%input is the meshgrid.

xratio = 10;   %hence always 20 cells.
yratio = 2; 

[X,Y]=meshgrid(1:size(im,2),1:size(im,1));  %cords of im

l = 1 : round(size(im,2)/xratio) : size(im,2);
w = 1 : round(size(im,1)/yratio) : size(im,1);
l=[l size(im,2)];
w=[w size(im,1)];
z=1;

for i=1:length(w)-1
    for j=1:length(l)-1  %across length first
        
        x = X(w(i):w(i+1),l(j):l(j+1));
        y = Y(w(i):w(i+1),l(j):l(j+1));
        imageCell(z).Location = [x(:) y(:)];
        imageCell(z).isText = true;
  
        z=z+1;
    end
end