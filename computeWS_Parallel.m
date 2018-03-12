function label = computeWS_Parallel(baseLine_id,label)

baseLine = label.traj_textLines{baseLine_id}; % most confident text line (the longest one)

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

%--------------------optimization for each WS line-----------------

for i=1:length(label.traj_WSLines)    
    
shortLine = label.traj_WSLines{i};

% fit a straight line through shortline.
% compute the offset for baseLine by min the distance.
[~,~,coeff] = linearFit(shortLine);

dis_old = label.CorredtedImage_Size(1);

for offset=0:0.1:label.CorredtedImage_Size(1)

    if(baseLine_id>=i)  % going in reverse direction than normal
        offset = -offset;
    end

    % new points
    P = baseLine + [xnormal*offset ynormal*offset];

    % distance of the points to the line along their normal direction
    dis = abs((P(:,1).*coeff(1) - P(:,2) + ones(length(P),1)*coeff(2)))/sqrt(coeff(1)^2 + 1);
    dis_new = mean(dis);

    % condition to stop
    if(dis_new > dis_old)
        break;                  % we have found the offset
    end

    dis_old = dis_new;
end

if(i==1)
    short = baseLine + (offset-4).*[xnormal ynormal];  % adding some buffer for the top and bottom WS line
elseif(i==length(label.traj_WSLines))
    short = baseLine + (offset+4).*[xnormal ynormal];
else
    short = baseLine + offset.*[xnormal ynormal];
end

% make the pts smooth
% [~,~,coeffs] =  cubicFit(short); 
% [~,~,coeffs] = quadraticFitNew(short);
[~,~,coeffs] = linearFitNew(short);

x_values = (1:label.CorredtedImage_Size(2));              % extrapolation
% pts = cubicFit_evaluate(x_values,coeffs);
% pts = quadraticFit_evaluate(x_values,coeffs);
pts = linearFit_evaluate(x_values,coeffs);
% pts = cropLines(label.mask,pts);

label.traj_WS_coeff{i} = coeffs;
label.traj_WSLines{i}=interparc(150,pts(:,1),pts(:,2),'linear');

end




end