function [ mse] = mseFunc(p)
%  Constrained objective function
%   Detailed explanation goes here
global imgMd; % Measured image with resolution model
global gaussFilter;
% global rampF;
% global winF;
img=createActImg2D( p );
imgBlur=imfilter(img,gaussFilter,'same');
%%%  Windowed ramp filter %%%%
% rImgBlur=imfilter(imgBlur,rampF,'same');
% rImgBlur(rImgBlur<0)=0;
% wrImgBlur=imfilter(rImgBlur,winF,'same');
% tmp=(wrImgBlur-imgMd).^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmp=(imgBlur-imgMd).^2;
mse=mean(tmp(:));
end

