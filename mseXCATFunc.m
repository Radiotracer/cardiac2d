function [ mse] = mseXCATFunc(p)
%  Constrained objective function
%  
global imgMd;
global gaussFilter;
img=createXCATImg2D( p );
imgBlur=imfilter(img,gaussFilter,'same');
tmp=(imgBlur-imgMd).^2;
% tmp=(img-imgMd).^2;
mse=mean(tmp(:));
end
