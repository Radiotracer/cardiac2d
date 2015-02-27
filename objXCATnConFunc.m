function [ val] = objXCATnConFunc(p)
%  mse function + shape constraint
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
global weight; % weight of the shape constraint;

nseg=floor((numel(p)-5)/2);
% compute the variance of radii and thicknesses
radiusP=zeros(1,nseg);
thicknessP=zeros(1,nseg);
for k=1:nseg
    radiusP(k)=p(6+2*(k-1));
    thicknessP(k)=p(7+2*(k-1));
end;
val=mseXCATFunc(p)+weight*(myvar(radiusP)+var(thicknessP));

end