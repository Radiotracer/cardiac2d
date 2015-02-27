function [ val] = objConFunc(p)
%  mse function + shape constraint
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
global weight; % weight of the shape constraint;
global nseg;

% compute the variance of radii and thicknesses
radiusP=zeros(1,nseg);
thicknessP=zeros(1,nseg);
for k=1:nseg
    radiusP(k)=p(5+3*(k-1));
    thicknessP(k)=p(6+3*(k-1));
end;
val=mseFunc(p)+weight*(var(radiusP)+var(thicknessP));

end
