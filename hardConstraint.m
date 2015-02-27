function [ c,ceq ] = hardConstraint(p)
%  Parameters: 
% (a) Center :p(1),p(2)
% (b) Blood pool activity; p(3)
% (c) Background activity; p(4)
% (d) Endocardium radii and corresponding myocaridum thicknesses(#=nRad):
%     [p(5),p(6)] --> [p(5+2*(nRad-1), p(6+2*(nRad-1)]
% (e) Activities of myocardium segments (#=nRad):
%  p(4+2*nRad+1)-->p(4+2*nRad+nRad)

global nRad;

eqIdx=1;

% #1: Uniform myocardial radii
for k=1:(nRad-1)
    ceq(eqIdx)=p(5+(k-1)*2)-p(5+k*2);eqIdx=eqIdx+1;
end

% #2: Uniform myocardium thickness
for k=1:(nRad-1)
    ceq(eqIdx)=p(6+(k-1)*2)-p(6+k*2);eqIdx=eqIdx+1;
end

% #3:Activity
c=[];



end

