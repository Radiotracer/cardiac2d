function [ out ] = myvar( r )
% Calcluate the new variance 
%   Detailed explanation goes here
r=r(:);
pr=[r;r(1)];
dpr=diff(pr);
out=sum(dpr.^2)/numel(r);
end

