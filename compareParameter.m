function [ differencePara ] = compareParameter( measuredPara, truePara )
% Compare the absolute difference between the measured parameters and the
% true parameters after rearranging the parameters into four different
% groups--center, activities, radii and thicknesses.
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
    tmpMeasured(5+(k-1))=measuredPara(7+3*(k-1));
    tmpMeasured(5+nseg+(k-1))=measuredPara(5+3*(k-1));
    tmpMeasured(5+2*nseg+(k-1))=measuredPara(6+3*(k-1));
    tmpTrue(5+(k-1))=truePara(7+3*(k-1));
    tmpTrue(5+nseg+(k-1))=truePara(5+3*(k-1));
    tmpTrue(5+2*nseg+(k-1))=truePara(6+3*(k-1));
end
 fprintf('%.1f,',tmpMeasured);

differencePara=tmpMeasured-tmpTrue;

end

