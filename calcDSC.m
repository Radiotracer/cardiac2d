function [dsc, dm] = calcDSC( p1, p2)
% This function calculates the dice similarity coefficient (DSC) of the
% left ventricle myocardium binary mask images created from parameters.
% (a) Center :p(1),p(2)
% (b) Blood pool activity; p(3)
% (c) Background activity; p(4)
% For each segment in myocardium (8 segments):
% The 1st segment
%(d)Central point radius on endocardium  p(5) 
%(f)Thickness p(6)
%(e)Myocardium activity p(7) 
% The qth segment: p(5+3*(q-1):7+3*(q-1))
[vol1,mask1] = calcVolOfMyocardium( p1 );
[vol2,mask2] = calcVolOfMyocardium( p2 );
intersection=mask1.*mask2;
if (vol1+vol2)~=0
    dsc=2*sum(intersection(:))/(vol1+vol2);
else
    dsc=-1;
end;
dm=abs(mask1-mask2);
end

