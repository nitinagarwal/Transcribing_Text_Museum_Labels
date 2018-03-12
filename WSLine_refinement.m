function label = WSLine_refinement(label)
% after WS line extrapolation, refine the WS lines so that they are wrap
% around the text.

img=rgb2gray(label.img_Correct);
last_xcord = label.CorredtedImage_Size(2);
x_values = (1:last_xcord);  
offset = 0:.01:1;
sz=1;

baseLine = label.traj_WSLines{1}; % all WSlines have same normal pts

% compute the normalized normal for each pt in baseLine
for i=2:size(baseLine,1)-1

        x = baseLine(i-1:i+1,1);
        y = baseLine(i-1:i+1,2);
        
        B=[x-x(2) y-y(2)]';         
        [U,S,~]=svd(B);

        ratio=S(1,1)/S(2,2);            
      
        xnormal(i) = U(1,2); 
        ynormal(i) = U(2,2);
end

% for the first and last point
xnormal(1)=xnormal(2);
ynormal(1)=ynormal(2);
xnormal(i+1)=xnormal(i);
ynormal(i+1)=ynormal(i);

xnormal=xnormal';
ynormal=ynormal';



% move WS lines first (attempt1) ; can move the textLines as well.
for i=1:length(label.traj_textLines)
    
    top = label.traj_WSLines{i};
    bttm  = label.traj_WSLines{i+1};
    
    normal = bttm - top; % direction is top -> bttm

    last_text = label.traj_textLines{i}(end,:);  
    text_length = find(top(:,1) < last_text(:,1));
    
    dis = [];
    Up_pts = [];
    Dwn_pts = [];
    %for each top and bttm compute the new location
    for j=1:text_length(end)
          
        intensity = [];
        for k =1:length(offset)
            newpt = round(top(j,:) + offset(k).*normal(j,:));
            int = img(newpt(2),newpt(1));
            intensity= [intensity; int];
        end
            
%         figure; hold on
%         title(sprintf('line %d',j))
%         plot(intensity,'g*-','MarkerSize',1);

        %smoothening the intensity.
        intensity=[linspace(1,length(intensity),length(intensity))' intensity];
        newintensity = smoothlines(double(intensity));

%       plot(newintensity(:,1),newintensity(:,2),'k*-','MarkerSize',1);

        textPts = line_refinement(newintensity(:,2));
   
        if (~isnan(textPts))
%             plot(newintensity(textPts,1),newintensity(textPts,2),'m*','MarkerSize',4);
            Up_pts = [Up_pts; top(j,:) + offset(textPts(1)).*normal(j,:)];
            Dwn_pts = [Dwn_pts; top(j,:) + offset(textPts(2)).*normal(j,:)];
            dis = [dis; [norm(offset(textPts(1)).*normal(j,:)) norm(offset(textPts(2)).*normal(j,:))]];
        end
    end

    %clean up of pts
    Up_pts=removeOutliers(Up_pts);
    Dwn_pts=removeOutliers(Dwn_pts);
    
    % moving the other WS points (extrapolation)
    up_pts = top + median(dis(:,1)).*[xnormal ynormal];
    up_pts = cropLines(label.mask,up_pts);
    
    x_min = min(Up_pts(:,1));
    x_max = max(Up_pts(:,1));
    prev = up_pts((up_pts(:,1)<x_min),:);
    next = up_pts((up_pts(:,1)>x_max),:);
    
    Up_pts = [prev(1,:) ;prev(2:5:end,:); Up_pts; next(1:5:end-1,:); next(end,:)];
%     [Up_pts,~,Up_coeff] =  cubicFit(Up_pts);
    [Up_pts,~,Up_coeff] = quadraticFitNew(Up_pts);
    
    label.ref_traj_WS_coeff{sz} = Up_coeff;
    label.ref_traj_WSLines{sz}=interparc(150,Up_pts(:,1),Up_pts(:,2),'linear');
    sz=sz+1;
    
    dwn_pts = top + median(dis(:,2)).*[xnormal ynormal];
    dwn_pts = cropLines(label.mask,dwn_pts);
    
    x_min = min(Dwn_pts(:,1));
    x_max = max(Dwn_pts(:,1));
    prev = dwn_pts((dwn_pts(:,1)<x_min),:);
    next = dwn_pts((dwn_pts(:,1)>x_max),:);
    
    Dwn_pts = [prev(1,:) ;prev(2:5:end,:); Dwn_pts; next(1:5:end-1,:); next(end,:)];
%     [Dwn_pts,~,Dwn_coeff] =  cubicFit(Dwn_pts);
    [Dwn_pts,~,Dwn_coeff] = quadraticFitNew(Dwn_pts);
    
    label.ref_traj_WS_coeff{sz} = Dwn_coeff;
    label.ref_traj_WSLines{sz}=interparc(150,Dwn_pts(:,1),Dwn_pts(:,2),'linear');
    sz=sz+1;
    
    


    
    

    
    %smooth lines
%     [~,~,Up_coeff] =  cubicFit(Up_pts);
%     [~,~,Dwn_coeff] =  cubicFit(Dwn_pts);
%     
%     x_values = (1:last_xcord);                  % extrapolation
%     Up_pts = cubicFit_evaluate(x_values,Up_coeff);
%     Up_pts = cropLines(label.mask,Up_pts);
%     
%     Dwn_pts = cubicFit_evaluate(x_values,Dwn_coeff);
%     Dwn_pts = cropLines(label.mask,Dwn_pts);
    
%     [Up_pts,Up_coeff] =  smoothlines_Tianetal(Up_pts); 
%     [Dwn_pts,Dwn_coeff] =  smoothlines_Tianetal(Dwn_pts);       
%    traj_textLines =  smoothlines_Tianetal_evaluation(coeffs,x_axis,traj_textLines(:,1));
    
%     [trajnew,bestModel] = ransac(Up_pts);
%     plot(trajnew{i}(:,1),trajnew{i}(:,2),'k*-','LineWidth',0.2,'MarkerSize',1);hold on
%     plot(trajnew{i}(:,1),polyval(bestModel,trajnew{i}(:,1)),'b-','LineWidth',1)
    
    
    
end


end



