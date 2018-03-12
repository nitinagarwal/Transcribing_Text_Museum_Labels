function  [E,W,true_dir] = edge_connectivity1(ptset,label)
% we compute the edge connectivity between candidate pts 
% returns edges and their weights
% weights = alpha * distance + beta*slope

debug = true;
numpts = size(ptset,1);

TR = delaunayTriangulation(ptset);

% if(debug)
% figure,triplot(TR);hold on
% plot(ptset(:,1),ptset(:,2),'r*','MarkerSize',2);
% set(gca,'YDir','reverse')
% end

% building a graph with weights.
E = edges(TR);  % all the edges

%computing the median of the slope of all edges and max dis
for i=1:length(E)
   edge(i,:) = [(ptset(E(i,2),1) - ptset(E(i,1),1)) (ptset(E(i,2),2) - ptset(E(i,1),2))];
   dis(i) = sqrt((ptset(E(i,2),2) - ptset(E(i,1),2))^2 + (ptset(E(i,2),1) - ptset(E(i,1),1))^2);   
end

%normalize edge vectors
N = arrayfun(@(n) norm(edge(n,:)), 1:size(edge,1));
N=N';
edge = [edge(:,1)./N edge(:,2)./N];

%true direction
long_edge = label.edges{4};     % bottom label edge
[npt,~,~]=linearFit(long_edge);
true_dir = (npt(1,:)-npt(end,:));
true_dir = true_dir./norm(true_dir);

% computing the weights of each edge
for i=1:length(edge)
    
    delta(i) = abs(dot(true_dir,edge(i,:)));

    if(dis(i)<1)
        W(i) = 0;
    else
         W(i) =  dis(i) * (1-delta(i));           % dis(i) *       %exp(slope(i)^1.5) ;
    end
    
end
W = W';
distance = dis;

% apart from delaunay tringulation add more edges . for every pt check a
% closest radius pts and if not already present add those edges

radius = 20; 
[N,D]=knnsearch(ptset,ptset,'k',radius);

for i=1:numpts
   
    edge = [ones(radius,1)*i N(i,:)'];
    dis  = D(i,:)';
    edge(1,:) = [];         % removing self edge
    dis(1) = [];
    
    id1 = ~ismember(edge,E,'rows');
    id2 = ~ismember(edge,[E(:,2) E(:,1)],'rows');
    id=and(id1,id2);
    edge = edge(id,:);      %edges not in delanuay
        
    % dis between those pts
    dis = dis(id); 

    new_edges  = [(ptset(edge(:,2),1) - ptset(edge(:,1),1)) (ptset(edge(:,2),2) - ptset(edge(:,1),2))];
    
    %normalize edge vectors
    nm = arrayfun(@(n) norm(new_edges(n,:)), 1:size(new_edges,1));
    nm=nm';
    new_edges = [new_edges(:,1)./nm new_edges(:,2)./nm];

    new_delta = arrayfun(@(n) abs(dot(true_dir,new_edges(n,:))), 1:size(new_edges,1));
    delta = [delta new_delta];
    distance = [distance dis'];
        
    weight =  dis .* (1-new_delta');  % dis .*
    
    weight(find(dis)<1)=0;      % do not trust very small dis
    
    W = [W ; weight];
    E = [E ; edge];
end

figure,subplot(3,1,1),histogram(delta,30),title('abs(cos(two vectors))'),...
    subplot(3,1,2),histogram(distance,'binWidth',2),title('distance'),...
    subplot(3,1,3),histogram(W,'binWidth',2),title('final weight')


end









