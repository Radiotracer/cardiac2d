function [dsc, dm] = calcDSC( p1, p2)
% This function calculates the dice similarity coefficient (DSC) of the
% left ventricle myocardium binary mask images created from parameters.
%  Parameters: 
% (a) Center :p(1),p(2)
% (b) Blood pool activity; p(3)
% (c) Background activity; p(4)
% (d) Endocardium radii and corresponding myocaridum thicknesses(#=nRad):
%     [p(5),p(6)] --> [p(5+2*(nRad-1), p(6+2*(nRad-1)]
% (e) Activities of myocardium segments (#=nSeg):
%  p(4+2*nRad+1)-->p(4+2*nRad+nSeg)

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

