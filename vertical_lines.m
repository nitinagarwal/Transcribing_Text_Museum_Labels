function [line1,line2]= vertical_lines(lines,image,flag)

warning off
% for i=1:length(lines)
%   theta(i)=lines(i).theta;
%   pt1(i,:)=lines(i).point1;
%   pt2(i,:)=lines(i).point2;
%   rho(i)=lines(i).rho;
% end
%  
% 1:center_value   vertical lines
% center_value+1:end horizontal lines

% computing the midpoint.
for i=1:length(lines)
    midpts(i,:)= [(lines(i).point2(1,1)-lines(i).point1(1,1))/2 ...
        (lines(i).point2(1,2)-lines(i).point1(1,2))/2];
    midpts(i,:)= [midpts(i,1)+lines(i).point1(1,1) midpts(i,2)+lines(i).point1(1,2)]; 
    
end

% figure,imshow(image)
% plot(midpts(:,1),midpts(:,2),'*m')

% Computing the farthest lines apart (with some deviation delta_distance)
array=midpts(:,:);                      %horizontal lines

combi=combntns(1:length(array),2);  % all combinations for 2 points. 

% computing distance between all the midpoints. 
for i=1:length(combi)
    dis(i)= sqrt((array(combi(i,1),1)-array(combi(i,2),1))^2 + ...
                (array(combi(i,1),2)-array(combi(i,2),2))^2);        
end

% computing the angle difference between all those lines. 
for i=1:length(combi)
    pt1=lines(combi(i,1)).point1;
    pt2=lines(combi(i,1)).point2;
    slope1 = abs(atand((pt2(1,2)-pt1(1,2))/(pt2(1,1)-pt1(1,1))));
    
    pt1=lines(combi(i,2)).point1;
    pt2=lines(combi(i,2)).point2;
    slope2 = abs(atand((pt2(1,2)-pt1(1,2))/(pt2(1,1)-pt1(1,1))));
    
    theta(i) = abs(slope1-slope2);
    
end
clear pt1 pt2

[val,~]=max(dis);
delta_distance=15;
delta_theta=25;

% computing all lines with some delta_distance error but having approx same angle
id=find(dis>=val-delta_distance & dis<=val+delta_distance & theta <= delta_theta);

% candidate midpoints.
pt1(:,1:2)=array(combi(id,1),:);
pt2(:,1:2)=array(combi(id,2),:);

% flipping pt1 and pt2 if the order is reverse. 
% lines
for i=1:size(pt1,1)
    if(pt1(i,1) >= pt2(i,1))     % switch in case of horizontal lines
            temp=pt1(i,:);
            pt1(i,:)=pt2(i,:);
            pt2(i,:)=temp;
    end
    
%     if(pt1(i,2) >= pt2(i,2))     % switch in case of vertical lines 
%             temp=pt1(i,:);
%             pt1(i,:)=pt2(i,:);
%             pt2(i,:)=temp;
%     end
    
end

if(flag)
    %plot all candidate distance lines.
    figure,imshow(image); hold on
    plot([pt1(:,1)';pt2(:,1)'],[pt1(:,2)';pt2(:,2)'],'-m','LineWidth',2);
end
    
% removing duplicate points if any.
pt1=unique(pt1,'rows');
pt2=unique(pt2,'rows');

%getting the longest length line in pt1 and pt2 (this is repetition of the
%pervious command)
[~,idx1]=ismember(pt1,midpts,'rows');
[~,idx2]=ismember(pt2,midpts,'rows');

% computing the length of all the lines from hough transform
for i=1:length(lines)
    distance(i) = sqrt( (lines(i).point1(1,1)-lines(i).point2(1,1))^2 + ...
                                (lines(i).point1(1,2)-lines(i).point2(1,2))^2 );
end

vec1=distance(idx1);
vec2=distance(idx2);

[~,idx1]=ismember(max(vec1),distance);
[~,idx2]=ismember(max(vec2),distance);

% points returns row wise
line1=[lines(idx1).point1(1,1) lines(idx1).point1(1,2);lines(idx1).point2(1,1) lines(idx1).point2(1,2)];
line2=[lines(idx2).point1(1,1) lines(idx2).point1(1,2);lines(idx2).point2(1,1) lines(idx2).point2(1,2)];

end





