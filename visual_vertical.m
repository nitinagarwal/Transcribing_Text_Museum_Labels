function visual_vertical(pts,dirs,flag)

nBin=180;
edges = linspace(0, 180, nBin + 1);
edgesBoundary = (edges(1:end-1) + edges(2:end)) / 2;

if(ndims(pts)==3)
    pts=squeeze(pts);
end



if(flag ==true) % veritcal computed from Tian method
    dirs = (64-dirs)+32;
    plot_pts = line(edgesBoundary(dirs)',pts);
else
    plot_pts = line(dirs',pts);
end
 
plot([plot_pts(:,1) plot_pts(:,3)]',[plot_pts(:,2) plot_pts(:,4)]','b-','LineWidth',2);
plot(pts(:,1),pts(:,2),'y*','MarkerSize',2);

end

function plot_pts = line(slope,pt)

slope = slope - 90; % because the tangent slope is perpendicular to gradient
% slope=180-slope;         % very imp somehow the slopes are all flipped
slope=tand(slope);

c=pt(:,2)-slope.*pt(:,1);

% solving for x for two diff y: pt(2)-10, pt(2)+10
for i=1:length(c)
    x(i,1)=(pt(i,2)-5-c(i))/slope(i);
    x(i,2)=(pt(i,2)+5-c(i))/slope(i);
end

plot_pts=[x(:,1) pt(:,2)-5 x(:,2) pt(:,2)+5];

% if the angle is quite small the line could be very long. 
distance = sqrt( (plot_pts(:,1)-plot_pts(:,3)).^2 + (plot_pts(:,2)-plot_pts(:,4)).^2);

% ideally distance should be close to 10. 
id  = find(distance > 20);

for i=1:length(id)
% solving for y for two diff x: pt(1)-5, pt(1)+5
 plot_pts(id(i),2) = slope(id(i)) * (pt(id(i),1)-5) + c(id(i));
 plot_pts(id(i),4) = slope(id(i)) * (pt(id(i),1)+5) + c(id(i));
 
 plot_pts(id(i),1) = (pt(id(i),1)-5);
 plot_pts(id(i),3) = (pt(id(i),1)+5);   
end
    
       
end
