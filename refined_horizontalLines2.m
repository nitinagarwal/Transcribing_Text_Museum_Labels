function [textLine,plot2] = refined_horizontalLines2(intensity,raw_intensity)
% input: intensity - smoothedout version of the raw_intensity.
% output: indices of the upperpoints & lowerpoints of textlines

for i=2:length(intensity)-1
   
    x1=[i-1 intensity(i-1)];
    x2=[i intensity(i)];
    x3=[i+1 intensity(i+1)];
    
    vec1 = [(x1(1)-x2(1)) (x1(2)-x2(2))];
    vec2 = [(x3(1)-x2(1)) (x3(2)-x2(2))];
    
    a = dot(vec1,vec2)/(norm(vec1)*norm(vec2)) ;
    angle(i)=acosd(a);
    
end

candi=find(angle < 150); % the real index. 130 150 

if(numel(candi) < 3)  % not enough candidate points found
    textLine = {};
    plot2 = 0;
    return 
end

% finding the peaks and valleys in noisy intensity plot
candi_inten = intensity(candi);
min_in = min(candi_inten);
max_in = max(candi_inten);
av_cand = (min_in+max_in)/2;

isValley=zeros(length(candi),1); % flag to consider is that candi is in reallity a line.
isPeak=zeros(length(candi),1); % flag to consider is that candi is in reallity a line.
flag=true; % searching for peaks

start = find(candi_inten==min(candi_inten));
% bool(start)=true;
isValley(start)=true;

for i=start:-1:2
   %peak
    if(candi_inten(i-1) > candi_inten(i) && candi_inten(i-1)>av_cand && flag)
       isPeak(i-1)=true;
       flag=false;
    %valley
    elseif(candi_inten(i-1) < candi_inten(i) && candi_inten(i-1)<av_cand && ~flag)
       
       isValley(i-1)=true;
       flag=true;      
    end
end
flag=true;
for i=start:1:length(candi)-1
   %peak
    if(candi_inten(i+1) > candi_inten(i) && candi_inten(i+1)>av_cand && flag)
       isPeak(i+1)=true;
       flag=false;
    %valley
    elseif(candi_inten(i+1) < candi_inten(i) && candi_inten(i+1)<av_cand && ~flag)
       
       isValley(i+1)=true;
       flag=true;      
    end
end

bool=xor(isPeak,isValley);
trajs_lines = candi((bool==1));
% trajs_lines = trajs_lines(2:2:end);  % just plotting the valleys for wacv

figure(50);hold on
plot2 = gca;
plot(plot2, intensity);hold on
plot(plot2, trajs_lines,intensity(trajs_lines),'g*');hold on
plot(plot2, raw_intensity,'r*-','MarkerSize',2)



if(isValley(1)==1) % valley cannot be on the first point
    isValley(1)=0;
    isPeak(2)=0;
end
if(isValley(end)==1) % valley cannot be on the last point
    isValley(end)=0;
    isPeak(end-1)=0;
end


% if any isValley==1 is left then proceed else return
if(numel(isValley==1)<1)
    textLine = {};
    plot2 = 0;
    return
end
   
% the valleys are the center of the textLines. Finding the textWidth.
% Assumption: the testWidth is equidistant from the center pt. (wrong
% assump)
id = find(isValley==1);
textCenter = candi(id);

textLine = [textCenter];







end












