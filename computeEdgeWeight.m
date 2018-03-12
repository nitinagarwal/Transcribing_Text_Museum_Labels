function cand_weig = computeEdgeWeight(cand_edg,ptset,true_dir,debug)

% debug=false;

% building a graph with weights.
E = cand_edg;  % all the edges

%computing the median of the slope of all edges and max dis
for i=1:length(E)
   edge(i,:) = [(ptset(E(i,2),1) - ptset(E(i,1),1)) (ptset(E(i,2),2) - ptset(E(i,1),2))];
   dis(i) = sqrt((ptset(E(i,2),2) - ptset(E(i,1),2))^2 + (ptset(E(i,2),1) - ptset(E(i,1),1))^2);
end

%normalize edge vectors
N = arrayfun(@(n) norm(edge(n,:)), 1:size(edge,1));
N=N';
edge = [edge(:,1)./N edge(:,2)./N];


for i=1:size(edge,1)
    
    delta(i) = abs(dot(true_dir,edge(i,:)));
    
    if(dis(i)<1)
        cand_weig(i)=0;
    else
        cand_weig(i) =  dis(i) * (1-delta(i)); 
    end
    
end
cand_weig = cand_weig';
   
if(debug)
figure,subplot(3,1,1),histogram(delta,30),title('abs(cos(two vectors))'),...
    subplot(3,1,2),histogram(dis,'binWidth',2),title('distance'),...
    subplot(3,1,3),histogram(cand_weig),title('final weight')
end



% 
% %computing the median of the slope of all edges and max dis
% for i=1:size(cand_edg,1)
%    sl(i,:) = [(ptset(cand_edg(i,2),2) - ptset(cand_edg(i,1),2))  (ptset(cand_edg(i,2),1) - ptset(cand_edg(i,1),1)) ];
%    dis(i) = sqrt((ptset(cand_edg(i,2),2) - ptset(cand_edg(i,1),2))^2 + (ptset(cand_edg(i,2),1) - ptset(cand_edg(i,1),1))^2);
% end
% sl = atan2d(sl(:,1),sl(:,2));
% sl((sl<0)) = sl((sl<0))+180;
% for i=1:length(sl)                      % computing the weights of each edge
%     slope(i) = abs(median_slope-sl(i));
%     if(slope(i)<10)
%     slope(i)=0;
%     end
%     cand_weig(i) = slope(i) * dis(i);        %dis(i) *          
% end

end
