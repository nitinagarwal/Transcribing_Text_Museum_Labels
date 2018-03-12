function [POINTS,RECTS] = compute_rect(gridUp,gridDwn,orien,sze)

pts = zeros(length(orien)*2,size(orien{1},2),2);  % total number of gridpts
pts(1,:,:) = gridUp{1};

%initial state
pt_location = squeeze(pts(1,:,:));
pt_orientation = orien{1};
line_segment = gridDwn{1};
pts(2,:,:) = inter(line_segment,pt_orientation,pt_location,sze);
newOrientation = interpolateOrientation(line_segment,orien{1},pts(2,:,:));
pt_orientation = newOrientation;

sz=3; % counter for keeping track for pts
for i=2:length(orien)
   
    pt_location = squeeze(pts(sz-1,:,:));
    line_segment = gridUp{i}; 
    pts(sz,:,:) = inter(line_segment,pt_orientation,pt_location,sze);
    
    newOrientation = interpolateOrientation(line_segment,orien{i},pts(sz,:,:));
    pt_orientation = newOrientation;
    sz=sz+1;
    
    pt_location = squeeze(pts(sz-1,:,:));
    line_segment = gridDwn{i};
    pts(sz,:,:) = inter(line_segment,pt_orientation,pt_location,sze);
    
    newOrientation = interpolateOrientation(line_segment,orien{i},pts(sz,:,:));
    pt_orientation = newOrientation;
    sz=sz+1;
    
end


% %initial state
% pt_location = squeeze(pts(1,:,:));
% pt_orientation = orien{1};
% % line_segment = gridUp{2};
% line_segment = gridDwn{1};
% 
% for i=2:length(orien)
%     
% pts(i,:,:) = inter(line_segment,pt_orientation,pt_location,sze); % passing the previous line, vertical ori, vertical pt loc     
% 
% newOrientation = interpolateOrientation(line_segment,orien{i},pts(i,:,:));
% 
% pt_location = squeeze(pts(i,:,:));
% pt_orientation = newOrientation;
% 
% if((i+1)<=length(gridUp))
%     line_segment = gridUp{i+1};
% else
%     break
% end
% 
% end

% % intersection on last line
% line_segment = gridDwn{end};
% pts(i+1,:,:) = inter(line_segment,pt_orientation,pt_location,sze);

POINTS = squeeze(pts(1,:,:));
for i=2:size(pts,1)
POINTS = [POINTS;squeeze(pts(i,:,:))];
end

% build rectangles
% order: left bottom -> right bottom -> right top -> left up..

num_row = size(pts,1);
num_col = size(pts,2);

RECTS = zeros( (num_row-1) * (num_col-1),4);  % #point * #WSlines
id=1;
for i=1:num_row-1
   for j=1:num_col-1
    
    RECTS(id,:) = [ i*(num_col)+j i*(num_col)+j+1 (i-1)*(num_col)+j+1 (i-1)*(num_col)+j];
    id=id+1;
   end
end

end

function newOrientation = interpolateOrientation(line,org_orien,query_pts)
% doing a linear interpolation of orientation
query_pts = squeeze(query_pts);

input = [line(:,1) org_orien'];

%fit a 15th order
coeffs=polyfit(input(:,1),input(:,2),4); %15th degree polynomial. 
% ynew=polyval(coeffs,input(:,1));
% newpoints = [input(:,1) ynew];

newOrientation=polyval(coeffs,query_pts(:,1));
newOrientation = round(newOrientation');

% newpoints = [query_pts(:,1) newOrientation]

end

function  intersections = inter(horizontal, orientation, location,sz)

% nBin=64;
nBin=180;
edges = linspace(0, 180, nBin + 1);
edgesBoundary = (edges(1:end-1) + edges(2:end)) / 2;

%modifying the horizontal so that all verticals can intersect
equ = line_equation(horizontal(1,:),horizontal(2,:));
first = [1 equ(1)+equ(3)];
equ = line_equation(horizontal(end-1,:),horizontal(end,:));
last = [sz(2) sz(2)*equ(1)+equ(3)];
horizontal = [first;horizontal;last];

hori_shited=circshift(horizontal,1,1);

for i=2:length(horizontal)                      % all the pieacewise lines formed from horizontal
lines(i,:) = line_equation(horizontal(i,:),hori_shited(i,:));
end
lines(1,:)=[];

for i=1:length(location)

% slope = edgesBoundary(orientation(i));
% slope=180-slope;                            % very imp somehow the slopes are all flipped
% slope=tand(slope);

slope = orientation(i) - 90; 
slope=tand(slope);

pt1 = [location(i,1) location(i,2)];
pt2 = [0 pt1(2)-40];
pt2(1) = pt1(1) + (pt2(2)-pt1(2))/slope;

l = line_equation(pt1,pt2);

[xpos,ypos] = compute (lines, l, horizontal);

if(xpos == Inf)
    error('Vertical orientation not intersection with previous line');
end

intersections(i,:) = [xpos ypos]; 

end

end

function [xpos,ypos] = compute(lines, l, horizontal)
% for a give vertical slope where does it intersect with set of lines

for i=1:length(lines)
    
    A = [-lines(i,1) -lines(i,2);-l(1) -l(2)];
    B = [lines(i,3); l(3)];
    X = linsolve(A,B);
    
    if ( X(1) >= horizontal(i,1) && X(1) < horizontal(i+1,1) ) % checking whether the intersection is true or not
        xpos = X(1);
        ypos = X(2);
        return;
    end
end

xpos=Inf;
ypos=Inf;

end

