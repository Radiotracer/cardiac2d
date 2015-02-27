function [img] = createActImg2D( p )
% This function creates a 256*256 2D image simulating a short-axis view of
% left ventricle (LV), modeled as a concentric circular sector with eight
% segments.
%  Parameters: 
% (a) Center :p(1),p(2)
% (b) Blood pool activity; p(3)
% (c) Background activity; p(4)
% For each segment in myocardium (8 segments):
% The 1st segment
%(d)Central point radius on endocardium  p(5) 
%(f)Thickness p(6)
%(e)Myocardium activity p(7) 
% The qth segment: p(5+3*(q-1):7+3*(q-1))
global dimX;
global dimY;
global nseg;
global dAng;
global hdAng;
inPts=zeros(2,nseg+1);
outPts=zeros(2,nseg+1);
% 
% [gx,gy]=meshgrid(1:dimX,1:dimY);
% ang=angle((gx-p(1))+1i*(gy-p(2)));
% nang=floor(ang/dAng)+5;   %% 1:nseg


[gx,gy]=meshgrid(1:dimX,1:dimY);
ang=angle((gx-p(1))+1i*(gy-p(2)));
ang(ang<0)=ang(ang<0)+2*pi;
nang=floor(ang/dAng)+1;   %% 1:nseg

for k=1:nseg
    ang=hdAng+dAng*(k-1);
    inPts(1,k)=p(1)+ p(5+3*(k-1))*cos(ang);
    inPts(2,k)=p(2)+ p(5+3*(k-1))*sin(ang);
    outPts(1,k)=p(1)+ (p(5+3*(k-1))+p(6+3*(k-1)))*cos(ang);
    outPts(2,k)=p(2)+ (p(5+3*(k-1))+p(6+3*(k-1)))*sin(ang);   
end
inPts(:,end)=inPts(:,1);
outPts(:,end)=outPts(:,1);
inCurve=fnplt(cscvn(inPts));
outCurve=fnplt(cscvn(outPts));


inMask=poly2mask(inCurve(1,:),inCurve(2,:),dimY,dimX);
outMask=poly2mask(outCurve(1,:),outCurve(2,:),dimY,dimX);
btwMask=outMask-inMask;
img=zeros(dimY,dimX);

for k=1:nseg
    tmp=nang;
    tmp(tmp~=k)=0;
    img=img+p(7+3*(k-1))/k*(btwMask.*tmp);
end;
img=img+p(3)*inMask+p(4)*(ones(dimY,dimX)-outMask);


%%% The activity assignment is not sufficiently efficient for 
%%% iterative search in optimization. Scanning pixels takes tremendous time, 
%%% especially the piecewise uniform myocardium activity. 
% for y=1:dimY
%     for x=1:dimX
%         if  ~inpolygon(x,y,outCurve(1,:),outCurve(2,:)) %Background
%             img(y,x)=p(4);                       
%         elseif inpolygon(x,y,inCurve(1,:),inCurve(2,:)) %Bloodpool
%             img(y,x)=p(3);
%         else  % Myocardium    
%             ang=angle((x-p(1))+1i*(y-p(2)));
%             if ang<0
%                 ang=ang+2*pi;
%             end;
%             angIdx=floor(ang/dAng);           
%             img(y,x)=p(7+angIdx*3);                
%         end;       
%     end;
% end

end

