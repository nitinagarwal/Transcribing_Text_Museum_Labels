function edgepoints = wavy_edgeDetection(edgeImage,point1,point2)
% computing the edges from two point set on the edgeImage

[M,N]=size(edgeImage);

B=1:M*N;
B=reshape(B,[N,M])';  % rowwise reshape

[A,xy]=img2graph_modified(B,edgeImage);
% [A,cost]=img2graph(B,edgeImage,false); hold on

% (columnindex1,rowindex1)
start_point = sub2ind([size(edgeImage,2),size(edgeImage,1)],point1(1),point1(2));
end_point = sub2ind([size(edgeImage,2),size(edgeImage,1)],point2(1),point2(2));

[PATH,FINALCOST] = dijkstra(A,start_point,end_point);

% [FINALCOST,PATH] = dijkstraold(A,xy,start_point,end_point);

[c,r]=ind2sub([size(edgeImage,2),size(edgeImage,1)],PATH);
% plot(c,r,'r*','MarkerSize',1);

% plot(c,r,'ks-','MarkerFaceColor','b'); hold on;
% plot(columnindex1,rowindex1,'ks-','MarkerFaceColor','y'); 
% plot(columnindex2,rowindex2,'ks-','MarkerFaceColor','y'); 

edgepoints=[c r];

end