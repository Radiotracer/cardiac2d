function [vol,mask] = calcVolOfMyocardium( p )
% This function calculates the volume of the myocardium in 
% a 256*256 2D image simulating a short-axis view of
% left ventricle (LV) created from parameter p.
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
vol=sum(btwMask(:));
mask=btwMask;


end

