%% Step 1: Computing the 4 corners and edges of the segmented image

% rectification of lytro images of labels.
% Note the input to the algorithm is a segmented label image.

close all;
clear all;

baseDir = pwd;  % Enter your base directory where images are stored
folderName = 'label/test_label'; %testing on test labels or real labels

dirinfo=dir(fullfile(baseDir,folderName,'/*_1.jpg'));           % reading  the images by name
% assert(length(dirinfo) == 3,' Less than 3 segmented labels ');

try
    load(fullfile(baseDir,folderName,'I'))
catch
    warning('did not load I. Might be running a new batch');
end


for image_num=1:length(dirinfo)

im = imread(fullfile(baseDir,folderName,dirinfo(image_num).name));    
    
tic    

label={};
label.img=im; 

% padding 
label.img = padarray(label.img,[50 50],255,'both');

label.debug=true;  % switch to generate images after each and every step. Useful for debuggig

% computing four corners
label = corner_detection(label);

% computing the edge pixels from the four corners
label = edge_detectionNew(label);

%% Step 2: Computing the pixels on the text

label = InterpolateLines3(label);

label.CorredtedImage_Size=[size(label.img_Correct,1) size(label.img_Correct,2)];

% subsampling the edges to 150 points
for i=1:4
    label.edges{i}=interparc(150,label.edges{i}(5:end-6,1),label.edges{i}(5:end-6,2),'spline');
end


% comuting the edge length of the label, whichever is max
left = abs(label.edges{2}(1,2) - label.edges{4}(1,2));
right = abs(label.edges{2}(end,2) - label.edges{4}(end,2));
label.width=round(max(left,right));

img=rgb2gray(label.img_Correct);

% show the edges
figure(1);
plot1=gca;
imshow(label.img_Correct,[]);hold on
color=['r';'g';'b';'m'];
for k=1:4
plot(label.edges{k}(:,1),label.edges{k}(:,2),[color(k) '*'],'MarkerSize',1);
end


candidatepts=[];

for i=1:length(label.trajAll_Vertical)
    
    pts=interparc(label.width,label.trajAll_Vertical{i}(5:end-5,1),label.trajAll_Vertical{i}(5:end-5,2)...
        ,'linear');
    
    intensity = [];
    
    for j=1:length(pts)  %removing the edge ptss
         int = img( round(pts(j,2)),round(pts(j,1)) );
         intensity= [intensity; int];  
    end

    h2=plot(plot1,pts(:,1),pts(:,2),'c*','MarkerSize',1);
    
    intensity=[linspace(1,length(intensity),length(intensity))' intensity];
    %smoothening the intensity.
    newintensity = smoothlines(double(intensity));
  
    [textLine,plot2] = refined_horizontalLines2(newintensity(:,2),intensity(:,2));
    
    if(~isempty(textLine))
        for g=1:length(textLine)
            plot(plot1, pts(textLine(g),1),pts(textLine(g),2),'r*','MarkerSize',2);
            candidatepts = [candidatepts; pts(textLine(g),:)];
        end
    end
     
    if(plot2~=0)
        cla(plot2);
    end

    delete(h2); %to clear the prev line
end


% removing false positive points
[gx,gy]=imgradientxy(img);
magSqr = sqrt(gx.^2 + gy.^2);

%removing the edges from magnitude. 
se = strel('disk',7);
mask=label.mask;
mask=imfill(mask);
mask = imerode(mask,se);

magSqr = magSqr.*mask;

pt = round(candidatepts);

for z=1:length(pt)
    in(z) = sum(magSqr(pt(z,2),pt(z,1)-2:pt(z,1)+2));
end
    
id = find(in<=300);
candidatepts(id,:)=[];
close all

figure,imshow(label.img_Correct,[]);hold on
plot(candidatepts(:,1),candidatepts(:,2),'r*','MarkerSize',1)



%% Step 3: Clustering the pts into textLines

% compute the edge connectivity between candidatepts.
[Edges,Weight,true_vec] = edge_connectivity1(candidatepts,label);
G = graph(Edges(:,1),Edges(:,2),Weight);

[T,pred] = minspantree(G);   % minimum spanning tree

[cand_edg(:,1),cand_edg(:,2)]=findedge(T);     % candidate edges to start optimization with
cand_weig = computeEdgeWeight(cand_edg,candidatepts,true_vec,true);

%visualize
figure(50); hold on
imshow(label.img_Correct,[]);hold on
plot(candidatepts(:,1),candidatepts(:,2),'r*','MarkerSize',2);hold on
set(gca,'YDir','reverse')
for i=1:length(cand_edg)
plot([candidatepts(cand_edg(i,1),1) candidatepts(cand_edg(i,2),1)], [candidatepts(cand_edg(i,1),2) candidatepts(cand_edg(i,2),2)],'g-','LineWidth',2);
end

%OPTIMIAZATION (remove edges with highest weight)

[cand_weig,id]=sort(cand_weig,'descend'); % edges sorted by edge weight
err = zeros(length(cand_weig),1);

bins = conncomp(T,'OutputForm','cell');   % for the initial graph
for j=1:length(bins)
    if(length(bins{j})>5)
       pts = candidatepts(bins{j},:); 
       [newpts,error,~] = cubicFit(pts);
        err(1)=err(1)+error;               % accumulating all errors
    end
end
plot(newpts(:,1),newpts(:,2),'k--','LineWidth',2);    

i=2;
while(i <= length(cand_weig) )  
% instead of optimizing remove all edges and compute error. 

Told{i}=T;  
%1. remove the edge with highest weight
    T = rmedge(T,id(1));          %always remove the top most weight
    
    %visualize
    plot([candidatepts(cand_edg(id(1),1),1) candidatepts(cand_edg(id(1),2),1)],...
        [candidatepts(cand_edg(id(1),1),2) candidatepts(cand_edg(id(1),2),2)],'m-','LineWidth',2);
    
%2. compute the connected components and fit a cubic spline to all comp
     bins = conncomp(T,'OutputForm','cell');
     
%3. compute the error (sum) for this fitting function 
    for j=1:length(bins)
        if(length(bins{j})>5)
           pts = candidatepts(bins{j},:); 
           [newpts,error,~] = cubicFit(pts);
            err(i)=err(i)+error;               % accumulating all errors
            ax1 = plot(newpts(:,1),newpts(:,2),'k--','LineWidth',2);
        end
    end
    
%4. stop if the error is not changing MUCH.
  if(err(i)==0)
      break;
  end
  delta(i) = err(i-1)/err(i) ;
  
  fprintf('iter %d..................has error %5.3f \n',i,err(i));
  disp(sprintf('delta is %4.3f \n',delta(i)));
  
% 5. Compute the edges and weight from fresh because edges in T are removed
    clear cand_edg cand_weig;
    [cand_edg(:,1),cand_edg(:,2)]=findedge(T);
    cand_weig = computeEdgeWeight(cand_edg,candidatepts,true_vec,false);
    
    [cand_weig,id]=sort(cand_weig,'descend');
    
i=i+1;
end
clear cand_edg cand_weig;

% to compute which iteration to pick, instead of picking the max(delta)
% pick max such that after that all values are very low
delta = delta(1:20);        

difference = abs(diff(delta)); % compute the derivative
idx=find(difference>0.3);      
idx=idx(end);

% there are more than 1 lines
Tf = Told{idx+1};

%visualize
figure,imshow(label.img_Correct,[]);hold on
plot(candidatepts(:,1),candidatepts(:,2),'r*','MarkerSize',1)
bins = conncomp(Tf,'OutputForm','cell');


% candidate pts are separated now
z=1;
for i=1:length(bins)
    if(length(bins{i})>5)         % removing those traj less than 5 pts 
        trajnew{z} = candidatepts(bins{i},:);
        plot(trajnew{z}(:,1),trajnew{z}(:,2),'m*-','LineWidth',0.2,'MarkerSize',1)
        z=z+1;
    end
end


%% ------------- RANSAC - remove outlier pts in each traj

for i=1:length(trajnew)
[trajnew{i},bestModel] = ransac(trajnew{i});
plot(trajnew{i}(:,1),trajnew{i}(:,2),'g*-','LineWidth',0.2,'MarkerSize',1);hold on
plot(trajnew{i}(:,1),polyval(bestModel,trajnew{i}(:,1)),'b-','LineWidth',1)
end


% sorting the textLines from top to bottom by y coordinate
for i=1:length(trajnew)
    ps(i) = trajnew{i}(1,2);
end
[~,id] = sort(ps);
temp = trajnew;
for i=1:length(trajnew)
trajnew{i} = temp{id(i)} ;
end

color=['r','g','b','c','y','m','k','w'];
figure,imshow(label.img_Correct,[]);hold on
for i=1:length(trajnew)
plot(trajnew{i}(:,1),trajnew{i}(:,2),color(i),'LineWidth',2)
end


Trajnew = trajnew;
% smoothening & resampling the textLine 150 points
for i=1:length(Trajnew)
    
    [traj_textLines,error,coeffs] =  cubicFit(Trajnew{i});       % fitting cubic function and not tial as ransac removes outliers
      
    label.traj_textLine_coeff{i}=coeffs;
    label.traj_textLinesOld{i} = traj_textLines;
    
    traj_textLines = unique(traj_textLines,'rows','stable');     % removes nan
 
    label.traj_textLines{i}=interparc(150,traj_textLines(:,1),traj_textLines(:,2),'linear');  % resampling the textLine
end

color=['r','g','b','c','y','m','k','w'];

%visualize
figure,imshow(label.img_Correct,[]);hold on
for i=1:length(Trajnew) 
plot(label.traj_textLines{i}(:,1),label.traj_textLines{i}(:,2),color(i),'LineWidth',2)
end


%% Step 4: Computing White Spaces

label = computeWSLines(label);

figure,imshow(label.img_Correct,[]);hold on
for i=1:length(label.traj_textLines) 
plot(label.traj_textLines{i}(:,1),label.traj_textLines{i}(:,2),'c-','LineWidth',1)
end

for i=1:length(label.traj_WSLines) 
plot(label.traj_WSLines{i}(:,1),label.traj_WSLines{i}(:,2),'r-','LineWidth',1)
end

close all;
%% Step 5: Extrapolating & Refining WSLines and textLines

for i=1:length(label.traj_textLines)
    x_cords(i) = label.traj_textLines{i}(end,1);
end

last_pt_id=find(x_cords == max(x_cords));


label = computeWS_Parallel(last_pt_id,label);

% visualization
figure,imshow(label.img_Correct,[]);hold on

for i=1:length(label.traj_WSLines) 
plot(label.traj_WSLines{i}(:,1),label.traj_WSLines{i}(:,2),'g-','MarkerSize',1)
end

% WSline refinement s.t. they enclose the text tightly
label = WSLine_refinement(label);

figure,imshow(label.img_Correct,[]);hold on
for i=1:length(label.ref_traj_WSLines)
plot(label.ref_traj_WSLines{i}(:,1),label.ref_traj_WSLines{i}(:,2),'g-','LineWidth',1)
end

% Recomputing the textLines from refined WSlines
label = textLine_refinement(label);

for i=1:length(label.traj_textLines) 
plot(label.ref_traj_textLines{i}(:,1),label.ref_traj_textLines{i}(:,2),'b-','LineWidth',1)
end

figure,imshow(label.img_Correct,[]);hold on
for i=1:length(label.ref_traj_textLines) 
   for j=1:length( label.ref_traj_textLines{i})
        if(label.pts_blank{i}(j)==0)    
            plot(label.ref_traj_textLines{i}(j,1),label.ref_traj_textLines{i}(j,2),'r*','MarkerSize',1)
        end
   end
end


%% Step 7: Vertical Orientation Estimation & Grid

label = vertical_line_9(label);  

flag=false;

figure,imshow(label.img_Correct,[]);hold on
for i=1:length(label.ref_traj_WSLines)
plot(label.ref_traj_WSLines{i}(:,1),label.ref_traj_WSLines{i}(:,2),'g-','LineWidth',1)
end

for i=1:length(label.verOri)
   visual_vertical(label.ver{i},label.verOri{i},flag);
end

label = build2Dgrid4(label);   % the orignal ver is better

figure,imshow(label.img_Correct,[]);hold on
vis_2d_rects(label.points, label.rect);

%% Step 8: Recification 

ver_num = label.gridpts_perLine;            % number of vertical lines
hor_num = length(label.ref_traj_WSLines);   % number of horizontal lines
pts = label.points;
rects = label.rect;

min_x = round(min(pts(:,1)));
max_x = round(max(pts(:,1)));
min_y = round(min(pts(:,2)));
max_y = round(max(pts(:,2)));

len = max_x-min_x;  % approx length and width of the label.
wid = max_y-min_y;

delta_x = round(len/ver_num);
delta_y = round(wid/hor_num);

xcords=reshape(pts(:,1),[ver_num,hor_num])';
ycords=reshape(pts(:,2),[ver_num,hor_num])';

% computing the aspect ratio for all rect
for i=1:hor_num-1
    for j=1:ver_num-1
        % compute aspect ratio in 2D grid
        l = sqrt((xcords(i,j) - xcords(i,j+1))^2 + (ycords(i,j) - ycords(i,j+1))^2) ;
        w = sqrt((xcords(i,j) - xcords(i+1,j))^2 + (ycords(i,j) - ycords(i+1,j))^2);
        ratio_rec(i,j) = l/w;   
    end
end

% matching coordinates in the final image with constraint = same aspect ratio as 2D grid
xcords_final = zeros(size(xcords));
ycords_final = zeros(size(ycords));

xcords_final(1:hor_num,1) = 1;   

for j=1:hor_num   % width of the horizontal lines would be approx same
   ycords_final(j,1:ver_num) = (j-1)*delta_y+1;
end

for i=1:hor_num-1
    for j=2:ver_num
        xcords_final(i,j) = xcords_final(i,j-1) + ratio_rec(i,j-1) * ...
            (ycords_final(i+1,j-1)-ycords_final(i,j-1)) ;               
    end
end

%averaging the xcords (because the xcords across column have to be same to
%keep the final a rectangle)
xcords_final = repmat(mean(xcords_final(1:hor_num-1,:),1),[hor_num 1]);
xcords_final = round(xcords_final);

[X,Y] = meshgrid(1:xcords_final(1,end),1:ycords_final(end,1));  % making the rectangle grid 
xq = zeros(size(X));  % these are the query pts whose values will be determined. 
yq = zeros(size(X));


for i=1:hor_num-1
    for j=1:ver_num-1       % comuting for each rect
        
       % top_left,top_right,btm_left, btm_right
       xq1 = [xcords(i,j) ycords(i,j)];         % orignal img
       xq2 = [xcords(i,j+1) ycords(i,j+1)];
       xq3 = [xcords(i+1,j) ycords(i+1,j)];
       xq4 = [xcords(i+1,j+1) ycords(i+1,j+1)];
       
       xq5 = [xcords_final(i,j) ycords_final(i,j)];   % new img
       xq6 = [xcords_final(i,j+1) ycords_final(i,j+1)];
       xq7 = [xcords_final(i+1,j) ycords_final(i+1,j)];
       xq8 = [xcords_final(i+1,j+1) ycords_final(i+1,j+1)];
       
       [x,y]=meshgrid(xcords_final(1,j):xcords_final(1,j+1),ycords_final(i,1):ycords_final(i+1,1));  % pts inside that rectangle
       
       transformed_pts= estimate_perspective([xq5;xq6;xq7;xq8],[xq1;xq2;xq3;xq4],[x(:) y(:)]);
       
       %reshape results
       xq(ycords_final(i,1):ycords_final(i+1,1),xcords_final(1,j):xcords_final(1,j+1)) = ...
           reshape(transformed_pts(:,1),size(x));
       yq(ycords_final(i,1):ycords_final(i+1,1),xcords_final(1,j):xcords_final(1,j+1)) = ...
           reshape(transformed_pts(:,2),size(y));
    end
end

[u,v] = meshgrid(-30:X(1,end)+30,-30:Y(end,1)+30); % for extrapolation

xq_new = interp2(X,Y,xq,u,v,'spline');
yq_new = interp2(X,Y,yq,u,v,'spline');

img = label.img_Correct;

clear recImg1;
recImg1(:,:,1) = interp2(im2double(img(:,:,1)),xq_new,yq_new);
recImg1(:,:,2) = interp2(im2double(img(:,:,2)),xq_new,yq_new);
recImg1(:,:,3) = interp2(im2double(img(:,:,3)),xq_new,yq_new);

recImg1 = imresize(uint8(recImg1*255),5);
close all;

label.rectPoints{1}=5*(xcords_final+30);  % padding of 30 pixels
label.rectPoints{2}=5*(ycords_final+30);
label.rectImage = recImg1;

figure,imshow(recImg1);
figure,imshow(img);

I{image_num} = label;
save(fullfile(baseDir,folderName,'I'),'I');  % saving all the variables

disp('processing done')
close all;

end


%% Step 9: Alignment of images from multiple views (run after all images are rectified)


val1 = 1;   % specify which images to align
val2 = 2;

label1=I{val1};
label2=I{val2};

I1=label1.rectImage;
I2=label2.rectImage;
figure,imshowpair(I1,I2,'falsecolor') % before alignment

%-----------------------affine alignment-----------------------

[I1, newimage,tform] = affineAlignment(I1,I2);

figure,imshowpair(I1,newimage,'falsecolor')

label1.affine = I1;
label2.affine = newimage;
label2.tform = tform;

% transforming the pts for the second image as well using tform
label2 = transformPoints(label2,tform);

%-----------------------non-linear alignment-----------------------
% non-linearly align each textLine so get accurate alignment

[label1,label2] = nonLinearAlignment(label1,label2);

I{val1}=label1;
I{val2}=label2;

save(fullfile(baseDir,folderName,'I'),'I');


%% Step 10: Compositing

a{1} = [];
a{2} = [];
%attempt 1 (taking max for each pixel)

try
a{1} = im2double(I{1}.aligned);
sz = [size(a{1},1) size(a{1},2)];
catch
end

try
a{2} = im2double(I{2}.aligned);
sz = [size(a{2},1) size(a{2},2)];
catch
end

A = zeros(sz(1),sz(2),3);
% taking max
if(~isempty(a{1}))
    A = bsxfun(@max,A,a{1});
end
if(~isempty(a{2}))
    A = bsxfun(@max,A,a{2}); 
end
figure,imshow(A); title('max value considered')

% attemp 2taking min
B = 255*ones(sz(1),sz(2),3);
% taking max
if(~isempty(a{1}))
    B = bsxfun(@min,B,a{1});
end
if(~isempty(a{2}))
     B = bsxfun(@min,B,a{2}); 
end
figure,imshow(B); title('min value considered')

% attempt 3 exposure fusion
ind=1;
if(~isempty(a{1}))
    combine(:,:,:,ind)=a{1};
    ind=ind+1;
end
if(~isempty(a{2}))
    combine(:,:,:,ind)=a{2};
    ind=ind+1;
end

C = exposure_fusion(combine,[1 1 0]);
figure,imshow(C);title('exposure fusion')

composition.max = A;
composition.min = B;
composition.expFus = C;

I{4} = composition;

save(fullfile(baseDir,folderName,'I'),'I');

disp('processing done')





