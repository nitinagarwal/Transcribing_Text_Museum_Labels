function textPts = line_refinement(intensity)
% input: intensity - smoothedout version of the raw_intensity.
% output: indices of the upperpoints & lowerpoints of textlines

%input intensity is all background. no text
if( var(intensity) < 1000 & min(intensity) > 150 )
    textPts = [nan nan];
    return;
end

for i=2:length(intensity)-1
   
    x1=[i-1 intensity(i-1)];
    x2=[i intensity(i)];
    x3=[i+1 intensity(i+1)];
    
    vec1 = [(x1(1)-x2(1)) (x1(2)-x2(2))];
    vec2 = [(x3(1)-x2(1)) (x3(2)-x2(2))];
    
    a = dot(vec1,vec2)/(norm(vec1)*norm(vec2)) ;
    angle(i)=acosd(a);
    
end

candi=find(angle < 170); % the real index. 150 130
candi(1)=[];             % this is always the first point
id = diff(candi);        % removing neighbouring pts (can come sometimes)
candi((id==1))=[];

if(numel(candi) < 3)  % not enough candidate points found
    textPts = [nan nan];
    return 
end



%----to ensure that there are two peaks at the end of the profile.

% candi has both the peaks and the valleys.
% finding the peaks and valleys in noisy intensity plot
isValley=zeros(length(candi),1); 
isPeak=zeros(length(candi),1); 

for i=2:length(candi)-1
   
    can = intensity(candi(i));
    prev = intensity(candi(i-1));
    after = intensity(candi(i+1));
    
    if( can > prev || can > after )
        isPeak(i)=true;
        isValley(i)=false;
        
    elseif(can < prev || can < after)
        isPeak(i)=false;
        isValley(i)=true;
    end
    
end

if(isPeak(2))       % for the boundary conditions
    isValley(1)=true;
elseif(~isPeak(2))
    isPeak(1)=true;
end

if(isPeak(end-1))
    isValley(end)=true;
elseif(~isPeak(end-1))
    isPeak(end)=true;
end

% ensuring there are peaks at both ends
if(isValley(1))
    candi = [1 candi];
    isValley = [0 ; isValley];
    isPeak = [1 ; isPeak];
end

if(isValley(end))
    candi = [candi length(intensity)];
    isValley = [isValley ; 0];
    isPeak = [isPeak ; 1];
end



% find the two points
candi_inten = intensity(candi);
min_in = min(candi_inten);
max_in = max(candi_inten);
av_cand = (min_in+max_in)/2;

diffs = candi_inten - av_cand;
id = find(diffs<0);

% an error condition
if(id(1)==1 || id(end) == length(candi))
    textPts = [nan nan];
    return;
end

%two peaks
first(1) = candi(id(1)-1);
second(1) = candi(id(end)+1);
%two valleys
first(2) = candi(id(1));
second(2) = candi(id(end));

textPts = round([mean(first) mean(second)]);  % taking the avg of first and second. kind of like non...
                                              % max suppression

end












