function label = build2Dgrid4(label)
% ver 4: interpolation of the orientation and computing intersection with
% above and below lines
% changed because textlines and WSlines are refined

% given a line and a slope compute the intersection if any.
sz=1;
for i=1:2:length(label.ref_traj_WSLines)

gridUp{sz} = inter(label.ref_traj_WSLines{i},label.verOri{sz},label.ver{sz},label.CorredtedImage_Size); % passing the previous line, vertical ori, vertical pt loc     
gridDwn{sz} = inter(label.ref_traj_WSLines{i+1},label.verOri{sz},label.ver{sz},label.CorredtedImage_Size);         
sz=sz+1;
end        

[pts,rects]=compute_rect(gridUp,gridDwn,label.verOri,label.CorredtedImage_Size);

label.points = pts;
label.rect = rects;

% label.gridUp{i}=gridUp{i};
% label.gridDwn{i}=gridDwn{i};
end

function  intersections = inter(horizontal, orientation, location,sz)

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








