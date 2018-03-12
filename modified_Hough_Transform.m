function lines = modified_Hough_Transform(clean_edge_Image,flag)
% Input is a relatively clean edge image. 

% Computing Hough Transform
[H,T,R] = hough(clean_edge_Image,'RhoResolution',0.5,'ThetaResolution',0.5);

%number of points for each of the three ranges
center_value=10;
edge_value=5;
range=100;  % i have changed this value from 140.(+/- 70 to +/- 50)

% for vertical lines (-25d to +25d)
center = H(:,181-range:181+range);
sortedindex1 = houghpeaks(center,center_value,'threshold',0);
for i=1:center_value
sortedindex1(i,2)=sortedindex1(i,2)+181-range-1;
sortedvalues1(i) = H(sortedindex1(i,1),sortedindex1(i,2));
end


P=sortedindex1;
V=sortedvalues1;
V=V./min(V(:));

if(flag)
    figure,imshow(imadjust(mat2gray(H)),'XData',T,'YData',R,...
          'InitialMagnification','fit');
    xlabel('\theta'), ylabel('\rho');
    axis on, axis normal, hold on;
    colormap(hot);
    scatter(T(P(1:center_value,2)),R(P(1:center_value,1)),round(V(1:center_value)'*100),'sg')
end

lines = houghlines(clean_edge_Image,T,R,P);

end









