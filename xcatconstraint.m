function [ c,ceq ] = xcatconstraint(p)
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

% ceq=zeros(3,1);
% ceq(1)=p(3)-17;
% ceq(2)=p(4)-5;
% ceq(3)=p(5)-93;
% ceq(4)=p(1)-71;
% ceq(5)=p(2)-65;
ceq=[];

nseg=floor((numel(p)-5)/2);
c=zeros(1,numel(p)*2);

% 50<p(1),p(2)<90
c(1)=p(1)-90;
c(2)=50-p(1);
c(3)=p(2)-90;
c(4)=50-p(2);

% 10<p(3)<25
c(5)=10-p(3);
c(6)=p(3)-25;

% 0<p(4)<10
c(7)=p(4)-10;
c(8)=-p(4);

% 80<p(5)<100
c(9)=80-p(5);
c(10)=p(5)-100;



% 10< radius <30, 8<thickness<12
for k=1:nseg
    c(11+4*(k-1))=10-p(6+2*(k-1));
    c(12+4*(k-1))=p(6+2*(k-1))-30;
    c(13+4*(k-1))=5-p(7+2*(k-1));
    c(14+4*(k-1))=p(7+2*(k-1))-15;
end

end
