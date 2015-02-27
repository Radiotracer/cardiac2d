close all;
data_pts=[1 sqrt(2)/2  0  -sqrt(2)/2 -1 -sqrt(2)/2 0 sqrt(2)/2 1;
    0 sqrt(2)/2  1  sqrt(2)/2  0  -sqrt(2)/2 -1 -sqrt(2)/2 0];
figure;
curve_pts=fnplt(cscvn(data_pts));
hold on;
axis equal;
plot(data_pts(1,:),data_pts(2,:),'ro'); hold off;


%% compare the fitted circle with the truth
circle_pts=[curve_pts(1,:);zeros(1,size(curve_pts,2))];
for k=1:size(circle_pts,2)
    if k<=(size(circle_pts,2)+1)/2
        circle_pts(2,k)=sqrt(1-circle_pts(1,k)^2);
    else
        circle_pts(2,k)=-sqrt(1-circle_pts(1,k)^2);
    end
end;
figure;plot(circle_pts(1,:),circle_pts(2,:),'ro'); axis equal;
figure;plot(curve_pts(1,:),curve_pts(2,:),'ro'); axis equal;
 mean(circle_pts(2,:)-curve_pts(2,:))
 
 %%
 close all;
 p=[128.5 128.5 20 5  50 30 80 50 30 85 50 30 90 50 30 95 50 30 100 ...
     50 30 105 50 30 110 50 30 115];
 
 nSeg=floor((length(p)-4)/3);
dAng=2*pi/nSeg;
inPts=zeros(2,nSeg+1);
outPts=zeros(2,nSeg+1);
for k=1:nSeg
    inPts(1,k)=p(1)+ p(5+3*(k-1))*cos(dAng*(k-1));
    inPts(2,k)=p(2)+ p(5+3*(k-1))*sin(dAng*(k-1));
    outPts(1,k)=p(1)+ (p(5+3*(k-1))+p(6+3*(k-1)))*cos(dAng*(k-1));
    outPts(2,k)=p(2)+ (p(5+3*(k-1))+p(6+3*(k-1)))*sin(dAng*(k-1));
end
inPts(:,end)=inPts(:,1);
outPts(:,end)=outPts(:,1);
fnplt(cscvn(inPts))
hold on;
axis equal;
fnplt(cscvn(outPts))

inCurve=fnplt(cscvn(inPts));
outCurve=fnplt(cscvn(outPts));

dimX=256;
dimY=256;
img=zeros(dimY,dimX);
for y=1:dimY
    for x=1:dimX
        if  ~inpolygon(x,y,outCurve(1,:),outCurve(2,:)) %Background
            img(y,x)=p(4);                       
        elseif inpolygon(x,y,inCurve(1,:),inCurve(2,:)) %Bloodpool
            img(y,x)=p(3);
        else  % Myocardium    
            ang=angle((x-p(1))+1i*(y-p(2)));
            if ang<0
                ang=ang+2*pi;
            end;
            angIdx=floor(ang/dAng);           
            img(y,x)=p(7+angIdx*3);                
        end;       
    end;
end
figure;imshow(img,[]);

%% Avoid loops; uniform myocardium activity
dimX=256;
dimY=256;
iBp=20;
iMy=80;
iBk=5;
inMask=poly2mask(inCurve(1,:),inCurve(2,:),dimY,dimX);
outMask=poly2mask(outCurve(1,:),outCurve(2,:),dimY,dimX);
img=(iBp-iMy)*inMask+(iMy-iBk)*outMask+iBk*ones(dimY,dimX);
figure;imshow(img,[]);

%% Avoid loops; nonuniform myocardium acitivity
close all; 
dimX=256;
dimY=256;
iBp=20;
iBk=5;
iMy=50*(1:nSeg);
[gx,gy]=meshgrid(1:dimX,1:dimY);
ang=angle((gx-p(1))+1i*(gy-p(2)));
figure;imshow(ang,[]);
nang=floor(ang/dAng)+nSeg;   %% nSeg is to avoid zero values
val=unique(nang(:)); 
figure;imshow(nang,[]);

inMask=poly2mask(inCurve(1,:),inCurve(2,:),dimY,dimX);
outMask=poly2mask(outCurve(1,:),outCurve(2,:),dimY,dimX);
btwMask=outMask-inMask;
img=zeros(dimY,dimX);

for k=1:nSeg
    tmp=nang;
    tmp(tmp~=val(k))=0;
    img=img+iMy(k)/val(k)*(btwMask.*tmp);
end;
img=img+iBp*inMask+iBk*(ones(dimY,dimX)-outMask);
figure;imshow(img,[]);



%% Debug the angle computation
close all;
dimX=256;
dimY=256;
cx=128;
cy=128;
nSeg=8;
dang=2*pi/nSeg;
hdang=pi/nSeg;
[gx,gy]=meshgrid(1:dimX,1:dimY);
ang=angle((gx-cx)+1i*(gy-cy));
ang(ang<0)=ang(ang<0)+2*pi;
figure;imshow(ang,[]);
nang=floor(ang/dang)+1;   %% nSeg is to avoid zero values
val=unique(nang(:))
numel(val)
figure;imshow(nang,[]);

%% Add noise & 2D ramp filter


close all;
tImg=imfilter(img,gaussFilter,'same');
figure;imshow(tImg,[]);
ntImg= double(1e+12*imnoise(tImg*1e-12,'poisson'));
figure;imshow(ntImg,[]);
[f1,f2] = freqspace(21,'meshgrid');
Hd = sqrt(f1.^2 + f2.^2);
% colormap(jet(64));
% mesh(f1,f2,Hd);
h = fsamp2(Hd);
% figure;freqz2(h);
fntImg=imfilter(ntImg,h,'same');
figure;imshow(fntImg,[]);

% tgaussFilter= fspecial('gaussian', [29 29],2);
% ffntImg=imfilter(fntImg,tgaussFilter,'same');
% figure;imshow(ffntImg,[]);

fntImg(fntImg<0)=0;
tgaussFilter= fspecial('gaussian', [29 29],2);
ffntImg_nz=imfilter(fntImg,tgaussFilter,'same');
figure;imshow(ffntImg_nz,[]);

figure;plot(tImg(:,128),'b-');
figure;plot(ntImg(:,128),'b-');
figure;plot(fntImg(:,128),'b-');
figure;plot(ffntImg_nz(:,128),'b-');


ftImg=imfilter(tImg,h,'same');
figure;imshow(ftImg,[]);
ftImg(ftImg<0)=0;
fftImg=imfilter(ftImg,tgaussFilter,'same');
figure;imshow(fftImg,[]);
tt=(ffntImg_nz-tImg).^2;
mean(tt(:))




close all;
tImg=imfilter(img,gaussFilter,'same');
figure;imshow(tImg,[]);
ntImg= double(1e+12*imnoise(tImg*1e-12,'poisson'));
figure;imshow(ntImg,[]);

[f1,f2] = freqspace(29,'meshgrid');
Hd = sqrt(f1.^2 + f2.^2);
colormap(jet(64));
mesh(f1,f2,Hd);
h = fsamp2(Hd);
figure;freqz2(h);
win= fspecial('gaussian', [29 29],3);
figure;freqz2(win);
winh=imfilter(h,win,'same');
figure;freqz2(winh);

fntImg=imfilter(ntImg,winh,'same');
figure;imshow(ntImg,[]);

%% Simulate noise: Add Gaussian noise to the 2D image and then convolve with the PSF
% SNR=1.4 u=93 ==> sigma=66.4
close all;
figure;imshow(img,[]);
nimg = imnoise(img,'gaussian');
nimg(nimg<0)=0;
figure;imshow(nimg,[]);
tnimg=imfilter(nimg,gaussFilter,'same');
figure;imshow(tnimg,[]);

%%
clear; close all;
gaussFilter= fspecial('gaussian', [29 29], 5.3);
tp=[128.5 128.5 17 5 30 40  80 30 40 93 30  40 93 30  40 93 30 40 93  30 40 93 30 40 93 30 40 93 ];
img=createActImg2D(tp);
figure;imshow(img,[]);

%%% Poisson noise %%%
nImg= double(1e+14*imnoise(img*1e-14,'poisson')); 
figure;imshow(nImg,[]);
tnImg=imfilter(nImg,gaussFilter,'same');
figure;imshow(tnImg,[]);

tImg=imfilter(img,gaussFilter,'same');
figure;imshow(tImg,[]);

tmp=(tImg-tnImg).^2;

%%
close all;
fid=fopen('lv50.raw');
imglv=fread(fid,[144 144],'float32','b')';
fclose(fid);
figure;imshow(imglv,[]);

%% 2/4/2015 
refM=double(imread('myo1.png'));
p=[74.0343   66.1840   17.7261    5.2462   85.3815   28.5217    7.3448   28.2145    8.3935   28.8888    8.1115   27.6556     8.1656   24.8689    8.4936   25.0295    8.9313   24.9838    8.2627   27.1252    9.8914   27.6679    9.1818   28.0081    7.9165   25.9052   10.6497   25.7741   10.3808];
[dsc, dm] = calcDSCX( p, refM);

p=[73.9342   67.1312   16.4918    4.9188   88.2494   27.6226    8.2472   27.7911    8.6142   28.1212    8.1212   26.6882    8.4398   24.3887    8.3625   24.5773    8.8909   24.6602    9.1063   27.6627    9.6872   28.8688    9.5697   28.1923    9.0283   26.1549   11.1056   25.9945   10.7902];
[dsc, dm] = calcDSCX( p, refM);

p=[74.0676   68.2913   17.3016    4.7762   87.3451   28.0668    7.4280   27.4469    7.2223   28.1662    6.4390   25.8143    9.5630   23.2933    7.2954   24.8771    9.0398   24.7104    8.6484   29.6967    8.4637   29.4543   10.0515   28.9092    9.8046   27.7890   10.1115   25.3209   11.3264];
[dsc, dm] = calcDSCX( p, refM);

p=[   70.9890   66.2507   14.3540    6.1793   90.3999   29.6416   10.1036   29.2913   12.9984   29.9305    9.1137   26.0139    7.9714   25.9395    5.6327   20.4294    9.5607   20.7023   11.3057   25.2416    8.2261   29.7655    7.2284   27.6130    8.9356   27.3574   7.9156   26.9581   13.7992];
[dsc, dm] = calcDSCX( p, refM);

p=[ 73.2707   66.5182   15.0728    6.1091   91.5360   27.6597    9.6155   27.6950   11.3660   29.7959    8.8816   27.3431    6.5488   25.9135    8.9063   24.4824    7.4568   23.8748    9.9460   26.9585    8.7046   29.1523    7.9222   28.0915    8.9241   26.3246  8.8191   25.7949   11.9190];
[dsc, dm] = calcDSCX( p, refM);

p=[74.0848   66.7860   14.8996    6.0970   91.1460   26.9317    9.7535   26.5281   11.7339   29.4361    8.7916   26.6203    7.8452   25.4809    9.6365   25.1907    7.9917   25.0650    9.2421   28.3281    7.6317   29.3722    7.8556   29.0261    8.4283   25.6610   8.1009   24.9259   11.6268];
[dsc, dm] = calcDSCX( p, refM);

p=[    72.0922   67.3545   14.7854    6.3312   92.0560   29.5579    9.0753   27.0618   12.4707   29.6106    8.5122   25.0256    8.5737   24.2727    9.0758   24.9004    5.9867   22.7792    9.1311   27.0323    8.2559   29.9376    7.2639   29.7323    8.9073   28.7986    7.7827   26.2434   13.5491];
[dsc, dm] = calcDSCX( p, refM);

p=[ 74.6579   66.8250   15.7955    6.2498   91.6070   27.1351    8.8673   27.0998   10.1511   28.6200    8.9978   27.4918  8.5631   26.6423    8.3529   26.4077    8.0965   26.4544    7.7075   28.5044    8.6482  29.1810    8.8520   28.0554  8.3901   25.4747    8.9716   25.7301   10.6508];
[dsc, dm] = calcDSCX( p, refM);

p=[74.5162   66.1837   15.6500    6.1794   90.5024   26.5159    9.9555   27.2575   10.3396   28.4383    9.9875   27.1285   8.9252   26.6206    9.1369   25.6171    8.9806   25.7321    9.1735   27.4872    9.2725   27.9459    9.3241   27.2895   9.2300   25.5320    8.9923   25.6073   10.2356];
[dsc, dm] = calcDSCX( p, refM);

p=[75.0373   66.8656   16.0429    6.2347   90.8429   26.6107    9.3994   26.9751    9.8853   27.4874   10.2085   26.7188    9.1996   26.5674    9.1786   26.0847    8.9454   26.3177    9.0768   27.0864    9.4781   27.5154    9.5944   26.9394    9.4289   25.8293    9.1461   26.1228    9.5584];
[dsc, dm] = calcDSCX( p, refM);

%% new constraint
p=[74.6824   65.9021   15.7360    6.1511   90.5724   26.8978    9.1560   28.1980    9.6462   28.8676    9.4647   27.9636    8.4424   26.9904    9.0605   26.4920    8.0117   26.6015    8.2756   27.3661    8.6675   27.8360    8.7032   27.1814    8.4873   25.8442   8.1041   25.7801    9.5410];
[dsc, dm] = calcDSCX( p, refM);

p=[74.8537   66.5357   15.6376    6.1897   89.7650   26.6615    9.3077   27.2993    9.8212   27.8111    9.9619   27.1544    9.1906   26.5040    9.2998   26.2365    8.7983   26.5389    9.1257   27.3950    9.4001   27.7451    9.5553   27.1720    9.2574   26.2229    8.7652   26.3580    9.6521];
[dsc, dm] = calcDSCX( p, refM);

p=[];
[dsc, dm] = calcDSCX( p, refM);

p=[];
[dsc, dm] = calcDSCX( p, refM);

p=[];
[dsc, dm] = calcDSCX( p, refM);

p=[];
[dsc, dm] = calcDSCX( p, refM);

%%


