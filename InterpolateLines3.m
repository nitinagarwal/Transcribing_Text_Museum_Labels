function label = InterpolateLines3(label)
% dense interpolation of lines along vertical direction
% Output: points inside the label

width = size(label.img_Correct,1);
len = size(label.img_Correct,2);

im = rgb2gray(label.img_Correct);
b=find(im~=255);  % points inside the label.


% vector of vertical edge
Vslope1 = [(label.cornerpts{1}(2,1)-label.cornerpts{1}(1,1)) (label.cornerpts{1}(2,2)-label.cornerpts{1}(1,2))];

Vslope = Vslope1; %abs(diff([Vslope1;Vslope2],1))/2;

traj={};
z=1;

% x-range for dense interpolation
[~,id] = min(label.cornerpts{1}(:,1));
pt1=label.cornerpts{1}(id,:);
[~,id] = max(label.cornerpts{2}(:,1));
pt2=label.cornerpts{2}(id,:);

s = (pt1(2) - 1) / Vslope(2);
beg = round(pt1(1) - s*Vslope(1));
s = (pt2(2) - 1) / Vslope(2);
en = round(pt2(1) - s*Vslope(1));


% connecting a line from top of the image to bottom of image with slope equal to label edge
for j=beg+5:2:en-5 % removing the edge cases

    step = (width - 0) / Vslope(2); % pt in top img (j,0)
    pt = [j 1;j+step*Vslope(1) width]; % pts in top and btm img 
    
    pts = round(interparc(round(1.5*width),pt(:,1),pt(:,2),'linear'));

    %need to check which pts lie inside the image
    id = find (pts(:,1) >= len | pts(:,1) <= 1 | pts(:,2) <= 1 | pts(:,2) >= width );
    if(~isempty(id))
        pts(id,:)=[]; % remove those points
    end
    
    % need to check which pts lie inside the label.
    can = sub2ind(size(im),pts(:,2),pts(:,1));
    an = intersect(can,b);
    
    if(~isempty(an))
        [r,c] = ind2sub(size(im),an);
        traj{z} = [c,r];
        z = z+1;
    end
    
end

%visualization    
if(label.debug)
    figure,imshow(label.img_Correct,[]);hold on
    for i=1:length(traj)
        plot(traj{i}(:,1),traj{i}(:,2),'*g','MarkerSize',1); 
    %      pause
    end
end

label.trajAll_Vertical = traj;

   
end
 