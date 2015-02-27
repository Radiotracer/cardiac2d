function [ val] = objConFunc(p)
%  mse function + shape constraint
%  Parameters: 
% (a) Center :p(1),p(2)
% (b) Blood pool activity; p(3)
% (c) Background activity; p(4)
% (d) Endocardium radii and corresponding myocaridum thicknesses(#=nRad):
%     [p(5),p(6)] --> [p(5+2*(nRad-1), p(6+2*(nRad-1)]
% (e) Activities of myocardium segments (#=nSeg):
%  p(5+2*nRad)-->p(4+nSeg+2*nRad)
global weight; % weight of the shape constraint;
global nRad;

% compute the variance of radii and thicknesses
radiusP=zeros(1,nRad);
thicknessP=zeros(1,nRad);
for k=1:nRad
    radiusP(k)=p(5+2*(k-1));
    thicknessP(k)=p(6+2*(k-1));
end;
val=mseFunc(p)+weight*(var(radiusP)+var(thicknessP));

end
