function [dsc, dm] = calcDSCX( p1, refM)
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

[vol1,mask1] = calcVolOfMyocardiumX( p1 );
intersection=mask1.*refM;
vol2=sum(refM(:));
if (vol1+vol2)~=0
    dsc=2*sum(intersection(:))/(vol1+vol2);
else
    dsc=-1;
end;
dm=abs(mask1-refM);
end