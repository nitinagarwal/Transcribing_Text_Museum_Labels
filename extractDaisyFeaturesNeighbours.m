function features = extractDaisyFeaturesNeighbours(cell2,id,p2,dzy2,flag)
% returns feture vector from the cells neighbouring to id. 
%if flag is true then neighbours are dwn else they are up.
% id ranges from 1:20

pts = [];
featureV = [];

iid = [id-1 id id+1];
if(flag)
    iid(iid<1 | iid>10)=[];        % its not cornercell
else
   iid(iid<11 | iid>20)=[];
end
    

for i=1:length(iid) 
    
    validpts = intersect(cell2(iid(i)).Location,p2,'rows');
    if(size(validpts,1)==0)
        continue;
    end
    featuresDaisy = extractDaisyFeatures(dzy2,validpts);
    pts = [pts;featuresDaisy.location];
    featureV = [featureV;featuresDaisy.vector];

end

if(flag)
    iid = [id+9 id+10 id+11];
    iid(iid<11 | iid>20)=[];        % its not cornercell
else
    iid = [id-9 id-10 id-11];
    iid(iid<1 | iid>10)=[];        % its not cornercell
end
    
for i=1:length(iid)
    
    validpts = intersect(cell2(iid(i)).Location,p2,'rows');
    if(size(validpts,1)==0)
        continue;
    end
    featuresDaisy = extractDaisyFeatures(dzy2,validpts);
    pts = [pts;featuresDaisy.location];
    featureV = [featureV;featuresDaisy.vector];

end

% remove duplicates if any
[pts,ia,~] = unique(pts,'rows'); 

features.vector = featureV(ia,:);
features.location = pts;

end
