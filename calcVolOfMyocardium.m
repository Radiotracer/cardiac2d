function [vol,mask] = calcVolOfMyocardium( p )
% This function calculates the volume of the myocardium in 
% a 256*256 2D image simulating a short-axis view of
% left ventricle (LV) created from parameter p.
%  Parameters: 
% (a) Center :p(1),p(2)
% (b) Blood pool activity; p(3)
% (c) Background activity; p(4)
% (d) Endocardium radii and corresponding myocaridum thicknesses(#=nRad):
%     [p(5),p(6)] --> [p(5+2*(nRad-1), p(6+2*(nRad-1)]
% (e) Activities of myocardium segments (#=nSeg):
%  p(4+2*nRad+1)-->p(4+2*nRad+nSeg)

global dimX;
global dimY;
global nRad;
global rAng; 
global hrAng; 

inPts=zeros(2,nRad+1);
outPts=zeros(2,nRad+1);
for k=1:nRad
    ang=hrAng+rAng*(k-1);
    inPts(1,k)=p(1)+ p(5+2*(k-1))*cos(ang);
    inPts(2,k)=p(2)+ p(5+2*(k-1))*sin(ang);
    outPts(1,k)=p(1)+ (p(5+2*(k-1))+p(6+2*(k-1)))*cos(ang);
    outPts(2,k)=p(2)+ (p(5+2*(k-1))+p(6+2*(k-1)))*sin(ang);   
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

