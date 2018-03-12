function [rotated_image,angle]= rotationCorrection1(image,flag)

mask = maskComputation(image);

[rows1,cols1]=find(mask==1);
r_max=max(rows1);
r_min=min(rows1);
c_max=max(cols1);
c_min=min(cols1);
ymid_image = r_min + (r_max-r_min)/2;
xmid_image = c_min + (c_max-c_min)/2;

B1=[cols1-xmid_image rows1-ymid_image]';            % centering the points
[U1,~,~]=svd(B1);                                   % SVD
image_angle=atand(U1(2,2)/U1(1,2));                 % orientation of the image

% removing some additional angle as it always underestimates the angles.
image_angle=image_angle-5;

if(image_angle > 0)
    angle=image_angle-90;
else 
    angle=90+image_angle;
end

new_image=padarray(image,[200 200],'replicate');          % padding instead of replicate making it 0 as some images are stuck to boundary    
new_image=imrotate(new_image,angle,'crop');
rotated_image=new_image(200:size(new_image,1)-200-1,200:size(new_image,2)-200-1,:);

if(flag)
    figure,
    subplot(1,2,1), imshow(image),hold on,
    plot(cols1,rows1,'c*','MarkerSize',1),
    title('Original Image');
    subplot(1,2,2), imshow(rotated_image),
    title('Rotated Image')
end

disp(['Rotation angle needed to Align is ' num2str(angle)])

end




