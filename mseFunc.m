function [ mse] = mseFunc(p)
%  Constrained objective function
%   Detailed explanation goes here
global imgMd; 
global gaussFilter;
img=createActImg2D( p );
imgBlur=imfilter(img,gaussFilter,'same');
tmp=(imgBlur-imgMd).^2;
mse=mean(tmp(:));
end

