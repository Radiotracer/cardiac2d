function [ c,ceq ] = hardConstraint(p)
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
% The qth segment: p(6+2*(q-1):7+2*(q-1))

global nseg;

eqIdx=1;

% #1: Uniform myocardial activity
for k=1:(nseg-1)
    ceq(eqIdx)=p(5+(k-1)*3)-p(5+k*3);eqIdx=eqIdx+1;
end

% #2: Uniform myocardium thickness
for k=1:(nseg-1)
    ceq(eqIdx)=p(6+(k-1)*3)-p(6+k*3);eqIdx=eqIdx+1;
end

% #3:Activity
c=[];



end

