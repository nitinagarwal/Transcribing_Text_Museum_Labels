function l = line_equation(pt1,pt2)
% % two points on a line and find the line equation
% pt1(x,y)
% pt2(x,y)
% equation of the form Ax

%Output of the form [slope -1 intercept];
%------------------ [xcoeff ycoeff intercept]  

if((pt2(1,1)-pt1(1,1)) == 0)   % vertical line
    
intercept = -pt2(1,1);
l=[1 0 intercept];

elseif((pt2(1,2)-pt1(1,2)) == 0)  % horizontal line

intercept= pt2(1,2);   
l=[0 -1 intercept];
    
else
    
slope = (pt2(1,2)-pt1(1,2))/(pt2(1,1)-pt1(1,1));
intercept = (pt1(1,2)*pt2(1,1)-pt2(1,2)*pt1(1,1))/(pt2(1,1)-pt1(1,1));

l=[slope -1 intercept];

end

end
