%% This study is to evaluate the resolution measurement error 
% True Parameter Settings
% Activity Ratio: Blood Pool : Myocardium : Background = 15 :95: 5
% Radius=20 mm, Thickness=8mm 
% For 256*256 matix size with 0.4mm*0.4mm pixel size,
% Radius=50, Thickness=20
% PET resolution FWHM=6mm, Sigma=6.0mm/(2sqrt(2ln2))=6.0/2.3548=2.5480mm=6.37p

%% Defect=50,  50% severity, 40 degree extent, weight=1;
tp_r=[128.5 128.5 15 5 ...
    50 20 50 20 50 20 50 20 50 20 50 20 50 20 50 20 ...
    95 95 95 95 95 95 95 50 50 95 95 95 95 95 95 95 95 95];

global imgMd;
global nRad; nRad=8;
global rAng; rAng=2*pi/nRad;
global hrAng; hrAng=pi/nRad;
global nSeg;nSeg=floor((length(tp_r)-4-2*nRad));
global sAng;sAng=2*pi/nSeg;
global dimX;dimX=256;
global dimY;dimY=256;
global weight;weight=1;
global gaussFilter;

img_r=createActImg2D(tp_r);
figure;imshow(img_r,[]);title('Truth(Resolution)');
trueGaussFilter=fspecial('gaussian', [29 29], 6.37);
initP_r=[127 126 20 8 ...
    45 27 46 26 44 28 45 30 47 26 46 25 43 29 44 30 ...
    70 70 70 70 70 70 70 70 70 70 70 70 70 70 70 70 70 70];

options_r = optimoptions(@fmincon,...
    'Display','iter',...
    'Algorithm','interior-point',...
    'FinDiffType','central',...
    'FinDiffRelStep',0.001,...
    'MaxFunEvals',10000 ...
    );

nNoise_r=40;
sigmas=[6.07 6.17 6.27 6.37 6.47 6.57 6.67];
nr=length(sigmas);
measuredImgs_r=zeros(dimX,dimY,nNoise_r);
pVals_r=zeros(nNoise_r,nr,length(tp_r));
fVals_r=zeros(nNoise_r,nr);
dscL_r=zeros(nNoise_r,nr);
%dmL_r=zeros(dimX,dimY,nr,nNoise_r);
dm=zeros(dimX,dimY);
for n=1:nNoise_r
    nImg_r= double(1e+15*imnoise(img_r*1e-15,'poisson')); 
    imgMd=imfilter(nImg_r,trueGaussFilter,'same');
    measuredImgs_r(:,:,n)=imgMd;
    for k=1:nr
        gaussFilter= fspecial('gaussian', [29 29], sigmas(k));
        [pVals_r(n,k,:),fVals_r(n,k)] = fmincon(@objConFunc,initP_r,[],[],[],[],[],[],@noconstraint,options_r); 
        [dscL_r(n,k), dm] = calcDSC(pVals_r(n,k,:), tp_r);
    end;  
end





% formats_r=['ro-';  'bo-'; 'go-'; 'mo-' ;'yo-'; 'co-';'r*-'];
% figure;plot(1:nSeg,tp_r(4+2*nRad+1:4+2*nRad+nSeg),'ko-');hold on;
% xlabel('Segment Index');ylabel('Activity Estimation');title('True Sigma=6.37');
% strLegend_r=cell(nr+1,1);
% strLegend_r{1}='truth';
% for k=1:nr
%     actL_r=pVals_r(k,4+2*nRad+1:4+2*nRad+nSeg);
%     plot(1:nSeg,actL_r,formats_r(k,:));
%     strLegend_r{k+1}=sprintf('Sigma=%.2f',sigmas(k));   
% end
% legend(strLegend_r);
% 
% hold off;
% figure; plot(sigmas,dscL_r,'b*-');
% xlabel('Measured Sigma');ylabel('DSC');title('Segmentation Results (Resolution)');
% 
% figure; plot(sigmas,fVals_r,'b*-');
% xlabel('Measured Sigma');ylabel('fVal');title('Final Objective Function(Resolution)');
