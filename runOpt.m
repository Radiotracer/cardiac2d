%%   True Parameter Settings
% Activity Ratio: Blood Pool : Myocardium : Background = 17 :93: 5
% Radius=12 mm, Thickness=16 mm 
% For 256*256 matix size with 0.4mm*0.4mm pixel size,
% Radius=30, Thickness=40
% PET resolution FWHM=5mm, Sigma=5mm/(2sqrt(2ln2))=5/2.3548=2.1233mm=5.3p
%%
clear; close all;
global gaussFilter;
gaussFilter= fspecial('gaussian', [29 29], 5.3);

% global gaussFilter1;
% gaussFilter1= fspecial('gaussian', [29 29], 5.67);

global imgMd;
tp=[128.5 128.5 17 5 30 40  80 30 40 93 30  40 93 30  40 93 30 40 93  30 40 93 30 40 93 30 40 93 ];
global nseg;nseg=floor((length(tp)-4)/3);
global dAng;dAng=2*pi/nseg;
global hdAng;hdAng=pi/nseg;
global dimX;dimX=256;
global dimY;dimY=256;
img=createActImg2D(tp);
figure;imshow(img,[]);
%%%%%% No noise %%%%
imgMd=imfilter(img,gaussFilter,'same');

%%% Poisson noise %%%
% tImg=imfilter(img,gaussFilter,'same');
% figure;imshow(tImg,[]);
% ntImg= double(1e+13*imnoise(tImg*1e-13,'poisson')); 
% figure;imshow(ntImg,[]);
% % [f1,f2] = freqspace(29,'meshgrid');
% % Hd = sqrt(f1.^2 + f2.^2);
% % global rampF;
% % rampF = fsamp2(Hd);
% % fntImg=imfilter(ntImg,rampF,'same');
% % fntImg(fntImg<0)=0;
% global winF;
% winF= fspecial('gaussian', [29 29],2); 
% imgMd=imfilter(ntImg,winF,'same');
% figure;imshow(imgMd,[]);

% nImg= double(1e+15*imnoise(img*1e-15,'poisson')); 
% figure;imshow(nImg,[]);
% imgMd=imfilter(nImg,gaussFilter,'same');
% figure;imshow(imgMd,[]);

%initP=[127 127 20 7 25 47 70 25 46 73 24  45 78 26  48 80 27 47 60  27 44 88 27 48 85 25 45 87 ];
%initP=[120 120 20 7 25 47 70 25 46 73 24  45 78 26  48 80 27 47 60  27 44 88 27 48 85 25 45 87 ];
initP=[120 120 20 7 26 53 70 25 56 73 22  49 78 31  48 80 29 53 60  30 44 88 19 55 85 25 45 87 ];
imgInit=createActImg2D(initP);figure;imshow(imgInit,[]);
%figure;imshow(imfilter(imgInit,gaussFilter,'same'),[]);
options = optimoptions(@fmincon,...
    'Display','iter',...
    'Algorithm','interior-point',...
    'FinDiffType','central',...
    'FinDiffRelStep',0.001,...
    'MaxFunEvals',10000 ...
    );
figure;
plot([mseFunc(initP) compareParameter(initP,tp)],'b*');
xlabel('Objective[1],Center[2-3],Activity[4-13],Radius[14-21],Thickness[22-29]');
ylabel('Difference');
title('Difference between Initial Parameters and Truth');


[ncep,ncfval] = fmincon(@mseFunc,initP,[],[],[],[],[],[],@noconstraint,options);
[ep,fval] = fmincon(@mseFunc,initP,[],[],[],[],[],[],@hardConstraint,options);


global weight;
weight=1;
[ep_w1,fval_w1] = fmincon(@objConFunc,initP,[],[],[],[],[],[],@noconstraint,options);

weight=0.1;
[ep_wd1,fval_wd1] = fmincon(@objConFunc,initP,[],[],[],[],[],[],@noconstraint,options);

weight=0.5;
[ep_wd5,fval_wd5] = fmincon(@objConFunc,initP,[],[],[],[],[],[],@noconstraint,options);

% calcDSC(tp,ep_wd1)
% calcDSC(tp,ep_wd5)
% calcDSC(tp,ep_w1)
% fval_wd1
% fval_wd5
% fval_w1


figure;
plot(0:30,zeros([1 31]),'k-');
hold on;
plot([mseFunc(ncep) compareParameter(ncep,tp)],'g+');
plot([mseFunc(ep_wd1) compareParameter(ep_wd1,tp)],'ro');
plot([mseFunc(ep_wd5) compareParameter(ep_wd5,tp)],'ko');
plot([mseFunc(ep_w1) compareParameter(ep_w1,tp)],'bo');
plot([mseFunc(ep) compareParameter(ep,tp)],'m*');
legend('Truth','Shape weight=0','Shape weight=0.1','Shape weight=0.5','Shape weight=1.0','Hard shape');
xlabel('Objective[1],Center[2-3],Activity[4-13],Radius[14-21],Thickness[22-29]');
ylabel('Difference with truth');

figure;
plot(0:30,zeros([1 31]),'k-');
hold on;
plot(compareParameter(ncep,tp),'g+');
plot( compareParameter(ep_wd1,tp),'ro');
plot( compareParameter(ep_wd5,tp),'ko');
plot(compareParameter(ep_w1,tp),'bo');
plot(compareParameter(ep,tp),'m*');
legend('Truth','Shape weight=0','Shape weight=0.1','Shape weight=0.5','Shape weight=1.0','Hard shape');
%xlabel('Center[1-2],Activity[3-12],Radius[13-20],Thickness[21-28]');
xlabel('C_x[1], C_y[2], A_{bp}[3], A_{bkg}[4],	A_{myoi}[5-12], R_i[13-20\], T_i[21-28]');
ylabel('Difference');
title('Difference:Estimated Parameters with Truth (Noise)');


ncepimg=createActImg2D(ncep);figure;imshow(ncepimg,[]);title('Shape weight=0');
ep_wd1_img=createActImg2D(ep_wd1);figure;imshow(ep_wd1_img,[]);title('Shape weight=0.1');
ep_wd5_img=createActImg2D(ep_wd5);figure;imshow(ep_wd5_img,[]);title('Shape weight=0.5');
ep_w1_img=createActImg2D(ep_w1);figure;imshow(ep_w1_img,[]);title('Shape weight=1');
epimg=createActImg2D(ep);figure;imshow(epimg,[]);title('Hard shape');
figure;imshow(img,[]);title('Truth');


figure;
plot(0:20,zeros([1 21]),'k-');
hold on;
plot(compareParameter1(ncep,tp),'r+');
plot( compareParameter1(ep_wd1,tp),'ko');
plot(compareParameter1(ep_w1,tp),'bo');
plot(compareParameter1(ep,tp),'m*');
legend('Truth','Shape weight=0.0','Shape weight=0.1','Shape weight=1.0','Hard shape');
%xlabel('Center[1-2],Activity[3-12],Radius[13-20],Thickness[21-28]');
xlabel('C_x[1], C_y[2],R_i[3-10], T_i[11-18]');
ylabel('Difference');

%title('Difference between Estimated Geometry Parameters and the Truth (Noisy)');

