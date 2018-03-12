function [dirs,diff] = interpolatingEdges(mag_hist,LDir,RDir,rank,nBin,blank,points)
% dirs = optimal orientation at that points
% diff = diff between the optimal and the desired

% computing everything in vector and then converting into degrees

LDir = [1 tand(LDir)];
LDir = LDir./norm(LDir);

RDir = [1 tand(RDir)];
RDir = RDir./norm(RDir);

edges = linspace(0, 180, nBin + 1);
edgesBoundary = (edges(1:end-1) + edges(2:end)) / 2;

mag_hist = squeeze(mag_hist);
points = squeeze(points);

% top rank values
[~,ind]=sort(mag_hist,1,'descend');

top_rank = ind(1:rank,:);
top_rank = edgesBoundary(top_rank);    % top 15 gradient orientations for all num_pts 

for i=1:length(points)
% for each pt, linearly interpolate the LDir and RDir 

t = (points(i,1) - points(1,1))/(points(end,1) - points(1,1));
desired_ori= LDir + t*(RDir-LDir);

if(blank(i)==0)
    
    desired_ori=atand(desired_ori(2)/desired_ori(1));  % converting back to degrees
    dirs(i) = desired_ori;
    diff(i) = 0;
    
else
    % find the closest among the top 15 values
    
    % convert to vector
    vec = [ones(rank,1) tand(top_rank(:,i))];
    nom_vec = bsxfun(@(A,B)(sqrt(A.^2+B.^2)),vec(:,1),vec(:,2));
    vec = [vec(:,1)./nom_vec vec(:,2)./nom_vec];
    
    for k=1:rank
        temp(k) = abs(dot(vec(k,:),desired_ori));
    end
    id = find(temp==max(temp)); 

    %     temp = abs(top_rank(:,i) - desired_ori);
    %     id = find(temp==min(temp));
    
    if(numel(id) > 1)
        id=id(1);
    end
    
    dirs(i) = top_rank(id,i);
    diff(i) = (1-temp(id));  % weights - lower weights mean accurate values
        
end
    
end




end