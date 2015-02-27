%% This study is to evaluate the severity of defect (0%,20%,40%,60%,80%,100%). 
% True Parameter Settings
% Activity Ratio: Blood Pool : Myocardium : Background = 20 :95: 5
% Radius=20 mm, Thickness=8mm 
% For 256*256 matix size with 0.4mm*0.4mm pixel size,
% Radius=50, Thickness=20
% PET resolution FWHM=5mm, Sigma=6.0mm/(2sqrt(2ln2))=6.0/2.3548=2.5480mm=6.37p

%% Defect=50,  50% severity, 40 degree extent
global gaussFilter;
gaussFilter= fspecial('gaussian', [29 29], 6.37);
global imgMd;
global nRad; nRad=8;
global rAng; rAng=2*pi/nRad;
global hrAng; hrAng=pi/nRad;
global dimX;dimX=256;
global dimY;dimY=256;
global nSeg;
global sAng;
global weight; weight=1.0;

tp_s=[128.5 128.5 15 5 ...
    50 20 50 20 50 20 50 20 50 20 50 20 50 20 50 20 ...
    95 95 95 95 95 95 95 95 95 95 95 95 95 95 95 95 95 95];
nSeg=floor((length(tp_s)-4-2*nRad));
sAng=2*pi/nSeg;

initP_s=[127 126 20 8 ...
    45 27 46 26 44 28 45 30 47 26 46 25 43 29 44 30 ...
    70 70 70 70 70 70 70 70 70 70 70 70 70 70 70 70 70 70];

options_s = optimoptions(@fmincon,...
    'Display','iter',...
    'Algorithm','interior-point',...
    'FinDiffType','central',...
    'FinDiffRelStep',0.001,...
    'MaxFunEvals',10000 ...
    );

trueGaussFilter=fspecial('gaussian', [29 29], 6.37);

nNoise_s=40;
ns=6;
severity=zeros(1,ns);
tps_s=repmat(tp_s,[ns 1]);
img_s=zeros(dimX,dimY,ns);
nImg_s=zeros(dimX,dimY);
pVals_s=zeros(nNoise_s,length(tp_s),ns);
fVals_s=zeros(ns,nNoise_s);
dscL_s=zeros(ns,nNoise_s);
%dmL_s=zeros(dimX,dimY,ns,nNoise_s);
dm=zeros(dimX,dimY);
measuredImgs_s=zeros(dimX,dimY,ns,nNoise_s);
for k=1:ns
    severity(k)=(k-1)*0.2;
    tps_s(k,28:29)=(tps_s(k,28:29)-tps_s(k,4))*(1-severity(k))+tps_s(k,4);
    img_s(:,:,k)=createActImg2D(tps_s(k,:));
    for n=1:nNoise_s
        nImg_s= double(1e+15*imnoise(img_s(:,:,k)*1e-15,'poisson'));
        imgMd=imfilter(nImg_s,trueGaussFilter,'same');    
        measuredImgs_s(:,:,k,n)=imgMd;
        [pVals_s(n,:,k),fVals_s(k)] = fmincon(@objConFunc,initP_s,[],[],[],[],[],[],@noconstraint,options_s); 
        [dscL_s(k), dm] = calcDSC(pVals_s(k,:), tps_s(k,:));       
    end
end

% formats_s=['ro-';  'bo-'; 'go-'; 'mo-' ;'yo-'; 'co-'];
% figure;plot(1:nSeg,tp_s(4+2*nRad+1:4+2*nRad+nSeg),'ko-');
% xlabel('Segment Index');ylabel('Activity Estimation');title('Severity Effects');hold on;
% strLegend_s=cell(ns+1,1);
% strLegend_s{1}='Baseline';
% for k=1:ns
%     actL_s=pVals_s(k,4+2*nRad+1:4+2*nRad+nSeg);
%     plot(1:nSeg,actL_s,formats_s(k,:));
%     strLegend_s{k+1}=sprintf('Severity=%.2f',severity(k));    
% end
% legend(strLegend_s);
% hold off;
% 
% figure; plot(severity,dscL_s,'b*-');
% xlabel('severity');ylabel('DSC');title('Segmentation Results');
% 
% figure; plot(severity,fVals_s,'b*-');
% xlabel('severity');ylabel('fVal');title('Final Objective Function(Severity)');