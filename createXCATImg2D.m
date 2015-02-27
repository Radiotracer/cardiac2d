function [img] = createXCATImg2D( p )
% This function creates a 144*144 2D image simulating a short-axis view of
% left ventricle (LV)
%  Parameters: 
% (a) Center :p(1),p(2)
% (b) Blood pool activity; p(3)
% (c) Background activity; p(4)
% (d) Myocardium activity: p(5)
% For each segment in myocardium (8 segments):
% The 1st segment
%(e)Central point radius on endocardium  p(6) 
%(f)Thickness p(7)
% The qth segment: p(6+2*(q-1):7+2*(q-1))

nseg=floor((numel(p)-5)/2);
dAng=2*pi/nseg;
hdAng=pi/nseg;
inPts=zeros(2,nseg+1);
outPts=zeros(2,nseg+1);
dimX=144;
dimY=144;

for k=1:nseg
    ang=hdAng+dAng*(k-1);
    inPts(1,k)=p(1)+ p(6+2*(k-1))*cos(ang);
    inPts(2,k)=p(2)+ p(6+2*(k-1))*sin(ang);
    outPts(1,k)=p(1)+ (p(6+2*(k-1))+p(7+2*(k-1)))*cos(ang);
    outPts(2,k)=p(2)+ (p(6+2*(k-1))+p(7+2*(k-1)))*sin(ang);   
end
inPts(:,end)=inPts(:,1);
outPts(:,end)=outPts(:,1);
inCurve=fnplt(cscvn(inPts));
outCurve=fnplt(cscvn(outPts));

inMask=poly2mask(inCurve(1,:),inCurve(2,:),dimY,dimX);
outMask=poly2mask(outCurve(1,:),outCurve(2,:),dimY,dimX);
btwMask=outMask-inMask;
img=p(3)*inMask+p(4)*(1-outMask)+p(5)*btwMask;
%figure;imshow(img,[]);

end

