function label = transformPoints(label,tform)
% transform pts using affine tform

T=tform.T';
sz = size(label.rectPoints{1});
x=label.rectPoints{1}(:);
y=label.rectPoints{2}(:);
z=ones(length(x),1);
oldcords = [x y z];
newcords = round(T*oldcords');
x=reshape(newcords(1,:),sz);
y=reshape(newcords(2,:),sz);
label.rectPoints{1} = x;
label.rectPoints{2} = y;


end