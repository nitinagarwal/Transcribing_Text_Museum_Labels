function modified_pts = cropLines(mask,input_pts)
% this function removes pts which are outside the label.

labelRegion = mask;
accessRegion=imfill(labelRegion);
id=1;

% modified_pts = [];
input_pts = abs(input_pts);

for i=1:length(input_pts)
    
    if(round(input_pts(i,2)) <1 || round(input_pts(i,1)) <1 || round(input_pts(i,2)) > size(labelRegion,1) ...
            || round(input_pts(i,1)) > size(labelRegion,2) )
        continue;
    end
    
    if(accessRegion(round(input_pts(i,2)), round(input_pts(i,1))))  % if inside the label
       modified_pts(id,:)=input_pts(i,:);
       id=id+1; 
    end
    
end



end