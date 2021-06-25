% single_site_test.m
%
% Author : C Fujiwara
%
% This is a very stupid piece of code which is meant to show how one goes
% from a raw image to a single site.

%% Files

% Raw data taken from the following folders
%file='Y:\Lattice Experiment Programs\Lattice Image Analysis\Lattice Test Images\B128x128\B_TestImage_Nph=165_FWHM=3.2px_128x128_0000.mat';
% file='Z:\Data\2016\2016 01 January 2016\25 January 2016\Spin_Mixture_Evaporated_Double_Gradient_Plane_Selection.mat';

f1='2018_PSF_Test.mat';
f2='Spin_Mixture_Evaporated_Double_Gradient_Plane_Selection.mat';

%% Load and Display raw Data
disp('loading data...');

data=load(f2);
images=data.images;

xL=[220 350];xL=[200 400];
yL=[190 300];yL=[150 350];
cL1=[0 300];
cL2=[0 50];

xL=[1 512];
yL=[1 512];
hf=figure;
clf
hf.Color='w';
hf.Name='Raw Data';
for kk=1:length(images)
   subplot(2,2,kk);
   Z=images{kk};
   imagesc(Z)
   axis equal tight
   caxis(cL1);
   xlim(xL);
   ylim(yL);
   colorbar
end
colormap(purplemap);

Z=images{1}+images{2};

Z=Z(yL(1):yL(2),xL(1):xL(2));
X=xL(1):xL(2);X=X';
Y=yL(1):yL(2);Y=Y';
%% point spread function
% psf='Y:\Lattice Experiment Programs\Lattice Image Analysis\GUI 1-2-2\Jan_29_PSF.mat';
psf=load(f1);
psf=psf.psf;

hf2=figure;
hf2.Name='Measured Point Spread Function';
clf
Zpsf=psf.psf;
Zpsf=Zpsf/max(max(Zpsf));
Xpsf=psf.x(1,:);
Ypsf=psf.y(:,1);

subplot(221);
imagesc(Xpsf,Ypsf,Zpsf)
axis equal tight
colormap(purplemap);
caxis([0 1]);

subplot(222);
plot(Xpsf,Zpsf(:,round(length(Ypsf)/2)));
colormap(purplemap);
caxis([0 1]);

subplot(223);
plot(Ypsf,Zpsf(round(length(Ypsf)/2),:));
colormap(purplemap);
caxis([0 1]);

% Get the sampling of the psf
dX=Xpsf(2)-Xpsf(1);
dY=Ypsf(2)-Ypsf(1);
dL=mean([dX dY]);

%% Image Parameters
mag=[82.6 83.2];         % Magnification x and y
pxsize=[16 16];          % Pixel size
pxsize_real=pxsize./mag; % um/pixel in real space
dPx=1e3*mean(pxsize_real);   % Size of pixel in nm
%% Rescale the image
sc=2;

Z=imresize(Z,sc);
X=imresize(X,[length(X)*sc 1]);
Y=imresize(Y,[length(Y)*sc 1]);

% After scaling, each pixel will represent a smaller distance
dPx=dPx/sc;
%% Gaussian PSF

% Gaussian radius of PSF in nm (Graham's Thesis);
sPSF=260;

% Choose the size of the PSF be 6*sigma
Npx=sPSF*6/dPx;
Npx=round(Npx);

% Make it an odd matrix so that zero is included.
Npx=(1-mod(Npx,2))+Npx;

fakepsf = fspecial('gaussian',Npx,sPSF/dPx);

xpsf=dPx*(1:Npx);
xpsf=xpsf-mean(xpsf); 


%% Deconvolution

% Define deconvolution
nLucy=5;
mypsf=fakepsf;

hf4=figure;
clf
hf4.Color='w';
hf4.Name='Deconvolved Data';
colormap(purplemap);
hf4.Position=[100 100 1200 600];

subplot(131);
imagesc(X,Y,Z);
axis equal tight
cc=[0 700];
cb=colorbar;
cb.Location='NorthOutside';
caxis(cc);
drawnow;
xlabel('px');
ylabel('px');

% Show the PSF
subplot(132)
cla
imagesc(xpsf,xpsf,mypsf);
axis equal tight
caxis([0 max(max(mypsf))]);
cb=colorbar;
cb.Location='NorthOutside';

drawnow;
xlabel('nm');
ylabel('nm');

% Perform the deconovlution
subplot(133)
out = richlucy(Z,mypsf,nLucy);
% out2=imgaussfilt(out,0.5*sc);
imagesc(X,Y,out);
axis equal tight

caxis([0 max(cc)*1.2]);
cb=colorbar;
cb.Location='NorthOutside';
disp('done');

