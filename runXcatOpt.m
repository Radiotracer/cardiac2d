%%   True Image Description
% Activity Ratio: Blood Pool : Myocardium : Background = 17 :93: 5
% For 144*144 matix size with 0.9766mm*0.9766mm pixel size,
% PET resolution FWHM=5mm, Sigma=5mm/(2sqrt(2ln2))=5/2.3548=2.1233mm=2.1742p
%%
clear; close all;
imgbkg=double(imread('bkg1.png'));
imgbp=double(imread('bp1.png'));
imgmyo=double(imread('myo1.png'));
img=imgmyo*93+imgbp*17+imgbkg*5;
%figure;imshow(img,[]);
global gaussFilter;
gaussFilter= fspecial('gaussian', [29 29], 2.2);
global imgMd;

nImg= double(1e+15*imnoise(img*1e-15,'poisson')); 
imgMd=imfilter(nImg,gaussFilter,'same');
figure;imshow(imgMd,[]);title('Measured Image');

initP=[71 65 17 5 93 25 10 25 10 25 10 25 10 25 10 25 10 25 10 25 10 25 10 25 10 25 10 25 10];
initImg = createXCATImg2D( initP );
%figure;imshow(initImg,[]);title('Initial Image');
options = optimoptions(@fmincon,...
    'Display','iter',...
    'Algorithm','interior-point',...
    'FinDiffType','central',...
    'FinDiffRelStep',0.01,...
    'MaxFunEvals',100000 ...
    );
[ncep,ncfval] = fmincon(@mseXCATFunc,initP,[],[],[],[],[],[],@xcatconstraint,options);
finalImg = createXCATImg2D( ncep );
figure;imshow(finalImg,[]);
figure;imshow(abs(finalImg-img),[]);

global weight;
weight=1;
[ep_w1,fval_w1] = fmincon(@objXCATConFunc,initP,[],[],[],[],[],[],@xcatconstraint,options);
Imgep_w1 = createXCATImg2D( ep_w1);
figure;imshow(Imgep_w1,[]);
figure;imshow(abs(Imgep_w1-img),[]);

weight=0.1;
[ep_wd1,fval_wd1] = fmincon(@objXCATConFunc,initP,[],[],[],[],[],[],@xcatconstraint,options);
Imgep_wd1 = createXCATImg2D( ep_wd1);
figure;imshow(Imgep_wd1,[]);
figure;imshow(abs(Imgep_wd1-img),[]);

weight=0.5;
[ep_wd5,fval_wd5] = fmincon(@objXCATConFunc,initP,[],[],[],[],[],[],@xcatconstraint,options);
Imgep_wd5 = createXCATImg2D( ep_wd5);
figure;imshow(Imgep_wd5,[]);
figure;imshow(abs(Imgep_wd5-img),[]);

weight=5;
[ep_w5,fval_w5] = fmincon(@objXCATConFunc,initP,[],[],[],[],[],[],@xcatconstraint,options);
Imgep_w5 = createXCATImg2D( ep_w5);
figure;imshow(Imgep_w5,[]);
figure;imshow(abs(Imgep_w5-img),[]);

weight=10;
[ep_w10,fval_w10] = fmincon(@objXCATConFunc,initP,[],[],[],[],[],[],@xcatconstraint,options);
Imgep_w10 = createXCATImg2D( ep_w10);
figure;imshow(Imgep_w10,[]);
figure;imshow(abs(Imgep_w10-img),[]);

weight=20;
[ep_w20,fval_w20] = fmincon(@objXCATConFunc,initP,[],[],[],[],[],[],@xcatconstraint,options);
Imgep_w20 = createXCATImg2D( ep_w20);
figure;imshow(Imgep_w20,[]);
figure;imshow(abs(Imgep_w20-img),[]);

%% New constraint 
weight=10;
[ep_nw10,fval_nw10] = fmincon(@objXCATnConFunc,initP,[],[],[],[],[],[],@xcatconstraint,options);
Imgep_nw10 = createXCATImg2D( ep_nw10);
figure;imshow(Imgep_nw10,[]);
figure;imshow(abs(Imgep_nw10-img),[]);

weight=20;
[ep_nw20,fval_nw20] = fmincon(@objXCATnConFunc,initP,[],[],[],[],[],[],@xcatconstraint,options);
Imgep_nw20 = createXCATImg2D( ep_nw20);
figure;imshow(Imgep_nw20,[]);
figure;imshow(abs(Imgep_nw20-img),[]);

