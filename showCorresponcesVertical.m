function showCorresponcesVertical(img1,img2,pt1,pt2,flag)
% combining the images and plotting point correspondences

assert(ndims(img1)==ndims(img2));

if(ndims(img1)==3)
    img1=rgb2gray(img1);
    img2=rgb2gray(img2);
end

if(size(img1,2) > size(img2,2))
    
    img2 = padarray(img2,[0 size(img1,2)-size(img2,2)],255,'post');

elseif(size(img1,2) < size(img2,2))
    
    img1 = padarray(img1,[0 size(img2,2)-size(img1,2)],255,'post');
        
end
   
bigImage = cat(1,img1,img2);
    
figure,imshow(bigImage);hold on

for i=1:length(pt1)
       x1plot(i)=pt1(i,1);
       x2plot(i)=pt2(i,1);
       y1plot(i)=pt1(i,2);
       y2plot(i)=pt2(i,2) + size(img1,1);
end

% i=1;
% plot([x1plot(i);x2plot(i)],[y1plot(i);y2plot(i)],'r-','LineWidth',1);hold on 
% plot(x1plot(i),y1plot(i),'go', 'MarkerSize',4);
% plot(x2plot(i),y2plot(i),'bo', 'MarkerSize',4);

% disp(sprintf('%i score is %2.4f',i,score(i)));

if(flag==true)
    
        [col,row]=ginput(2); 
        [val,id] = find( x1plot >=col(1) & x1plot <=col(2) & y1plot >=row(1) & y1plot <=row(2) );
    
        if (numel(val)~=0)
           
            plot([x1plot(id(1:end));x2plot(id(1:end))],[y1plot(id(1:end));y2plot(id(1:end))],'r-','LineWidth',1);hold on 
            plot(x1plot(id(1:end)),y1plot(id(1:end)),'go', 'MarkerSize',4);
            plot(x2plot(id(1:end)),y2plot(id(1:end)),'bo', 'MarkerSize',4);
            
        end
    
    
    
%     for i=2:length(x1plot)
%         
%         plot([x1plot(1:i-1);x2plot(1:i-1)],[y1plot(1:i-1);y2plot(1:i-1)],'y-','LineWidth',1);hold on
%         plot([x1plot(i);x2plot(i)],[y1plot(i);y2plot(i)],'r-','LineWidth',1);hold on 
%         plot(x1plot(i),y1plot(i),'go', 'MarkerSize',4);
%         plot(x2plot(i),y2plot(i),'bo', 'MarkerSize',4);
% %         disp(sprintf('%i score is %2.4f',i,score(i)));
%         pause
%     end
else

plot([x1plot;x2plot],[y1plot;y2plot],'y-','LineWidth',1);hold on 
plot(x1plot,y1plot,'go', 'MarkerSize',4);
plot(x2plot,y2plot,'bo', 'MarkerSize',4);

end



end