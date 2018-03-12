function label = textLine_refinement(label)

sz=1;

for j=1:2:length(label.ref_traj_WSLines)

    up  = label.ref_traj_WSLines{j};
    down = label.ref_traj_WSLines{j+1};
    
    
p = [up; down]; 
constraint(:,1:2) = [[1:(length(up)-1)]' [2:length(up)]'];
constraint=[constraint; [(length(up)+1:length(up)+length(down)-1)' (length(up)+2:length(up)+length(down))']];
constraint=[constraint; [1 (length(up)+1)];[length(up) (length(up)+length(down))] ];

TR = delaunayTriangulation(p,constraint);
tf = isInterior(TR);  % finding all the traingles which are inside the convex hull of the boundary
TR=triangulation(TR(tf, :),TR.Points); 

dis_up = sqrt( (up(1,1)-up(end,1))^2+(up(1,2)-up(end,2))^2);
dis_dwn = sqrt( (down(1,1)-down(end,1))^2+(down(1,2)-down(end,2))^2);

list = TR.ConnectivityList;
% removing the last triangle of the short line as the voronoi pts curl up
if ( dis_up < dis_dwn)
     
    pt = ismember(TR.Points,up(end,:),'rows');
    pt=find(pt==1); 
    id1=find(list(:,1)==pt);
    id2=find(list(:,2)==pt);
    id3=find(list(:,3)==pt);
    id = [id1;id2;id3];
    
else
    pt = ismember(TR.Points,down(end,:),'rows');
    pt=find(pt==1); 
    id1=find(list(:,1)==pt);
    id2=find(list(:,2)==pt);
    id3=find(list(:,3)==pt);
    id = [id1;id2;id3];
end

% sometimes with CDT, additional pts are added. these pts usually form
% triangle among themselves. trying to find these pts. 
p_prime = TR.Points;
id1 = ismember(p_prime,p,'rows');
id1 = find(id1==0);
% these extra pts are usually for the boundary
pts_look_up = [1:1:length(up)]';
pts_look_dwn = [length(up)+1:length(up)+length(down)]';

if (j==1)
    pts_look_up = [pts_look_up;id1];
elseif(j==length(label.traj_textLines)+1)
   pts_look_dwn = [pts_look_dwn;id1];
end


% %removing traingles formed from the same line
   id1 = ismember(list(:,1),pts_look_up);
   id2 = ismember(list(:,2),pts_look_up);
   id3 = ismember(list(:,3),pts_look_up);
   id4=id1+id2+id3;
   id5 = find(id4==3);
   id = [id;id5];
   
   id1 = ismember(list(:,1),pts_look_dwn);
   id2 = ismember(list(:,2),pts_look_dwn);
   id3 = ismember(list(:,3),pts_look_dwn);
   id4=id1+id2+id3;
   id5 = find(id4==3);
   id = [id;id5];
   

list(id,:)=[];
TR=triangulation(list,TR.Points);
        
%visualize
% triplot(TR);
% set (gca,'Ydir','reverse');


[CC,~] = circumcenter(TR); % circumcenters / voronoi vertices
ti = pointLocation(TR,CC); % removing all points corresponding to sliver triangles

z=~isnan(ti);                               % checking if the circumcenters are within ANY triagulation
CC(z==0,:)=[]; 

% plot(CC(:,1),CC(:,2),'g*','MarkerSize',1)

[~,order]=sort(CC(:,1)); % sorting x cord
CC = CC(order,:);

% [CC,coeffs] =  smoothlines_Tianetal(CC);
% [CC,~,coeffs] =  cubicFit(CC);

% resampling 
CC=interparc(150,CC(:,1),CC(:,2),'linear');
label.ref_traj_textLines{sz}=CC;


x_begin = label.traj_textLines{sz}(1,1);     % first pt on text line
x_end = label.traj_textLines{sz}(end,1);     % the last pt on text line

id_begin = find(CC(:,1)<(x_begin-1));
id_end = find(CC(:,1)>(x_end+1));

label.pts_blank{sz} = ones(150,1);
label.pts_blank{sz}(id_begin)=0;               % they are blank or extrapolated points
label.pts_blank{sz}(id_end)=0; 

sz=sz+1;


% label.traj_WS_coeff{j}=coeffs;
clear constraint;
clear id;

end



end





