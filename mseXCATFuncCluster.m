function [ mse] = mseXCATFuncCluster(p)
%  Constrained objective function
%  
global imgMd;
img=createXCATImg2D( p );

writeim('tmp.im',img);
[status1,cmdout1]=system('genprj genprj.par  tmp.im  tmpprj.im');
[status2,cmdout2]=system('osem osem.par tmpprj.im  tmpRecon');
imgBlur=readim('tmpRecon.5.im');

tmp=(imgBlur-imgMd).^2;
mse=mean(tmp(:));
end
