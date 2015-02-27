function [ differencePara ] = compareParameter1( measuredPara, truePara )
% Compare the difference between the measured parameters and the
% true parameters for geometrical parameters
% groups--center,radii and thicknesses.
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
global nseg;
tmpMeasured=measuredPara;
tmpTrue=truePara;
for k=1:nseg
    tmpMeasured(3+(k-1))=measuredPara(5+3*(k-1));
    tmpMeasured(3+nseg+(k-1))=measuredPara(6+3*(k-1));    
    tmpTrue(3+(k-1))=truePara(5+3*(k-1));
    tmpTrue(3+nseg+(k-1))=truePara(6+3*(k-1));
end
tmpMeasured=tmpMeasured(1:(2+2*nseg));
tmpTrue=tmpTrue(1:(2+2*nseg));
differencePara=tmpMeasured-tmpTrue;

end

