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

%% Load Raw Data 

data=load(f2);
images=data.images;
% 
% 
% f2='Y:\Data\2021\2021.06\06.24\iXonUltra_2021-06-24_21-46-02.mat';
% data=load(f2);
% data=data.data;
% images{1}=imgaussfilt(data.RawImages(:,:,2)-200,1);
% images{2}=imgaussfilt(data.RawImages(:,:,2)-200,1);
% images{3}=imgaussfilt(data.RawImages(:,:,1)-200,1);

xL=[220 350];xL=[1 512];
yL=[190 300];yL=[1 512];
cL1=[0 300];
cL2=[0 50];

hf=figure;
clf
hf.Color='w';
hf.Name='Raw Data';
for kk=1:length(images)
   subplot(2,2,kk);
   Z=images{kk};
   %Z=imgaussfilt(Z,.5);
   imagesc(Z)
   axis equal tight
   caxis(cL1);
   xlim(xL);
   ylim(yL);
   colorbar
end
colormap(purplemap);


% Z=images{1}+images{2};
% Z=imgaussfilt(Z,);

Z=images{2}-images{1};

%%
% Load point spread function
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
% caxis([0 300])
% %%
% hf3=figure(3);
% clf
% hf3.Color='w';
% hf3.Name='Crop and Subtracted Data';
% Z=images{1}-images{3};
% for kk=1:length(images)
%    %Z=imgaussfilt(Z,.5);
%    imagesc(Z)
%    axis equal tight
%    caxis(cL2);
%    xlim(xL);
%    ylim(yL);
% end
% colormap(purplemap);


%% PSF and Deconvoution
%
% The electric field from the atomic light emitters (effectively point
% sources) get mapped onto the image with the Point Spread Function.
% Because there are many single emitters, these means that the imaged
% intensity pattern is the convolution of the "true" light emitters pattern
% and the point spread functino.
%
% image = convolution(atoms,PSF)
%
% Via the convolution theorem, we then know that
%
% FT(image)=FT(atoms)*FT(PSF)
%
% Which suggests that the atomic signal can be obtained by
%
% atoms = FT_inv(FT(image)/FT(PSF))
%
% The fourier transform of the PSF is also known as the optical transfer
% fucntion (OTF). An ideal OTF is unity across all spatial frequencies, but
% any real system will have a OTF that dies off at higher frequencies.
%
% This is operation is fairly straightfoward with the exception of what to
% do when the OTF reaches zero or is dominated by noise.  This send the denominotrer of the FT space
% to zero.  Clearly this is non-physical. The majority of computations is
% basically finding the best way to deal with this fact.
%
%

hf5=figure;
clf
hf5.Color='w';
hf5.Name='PSF Example';

% Create a PSF
N=70; % size of PSF (NxN)
s=5;  % gaussian radius (pixels)
PSF = fspecial('gaussian',N,s);
PSF=PSF;

% Create an ideal image
I=zeros(200,200);
I(100,100)=150;
I(150,150)=150;
I(180,180)=150;
I(100,85)=150;

% Convolve the ideal image with the PSF
blurred = imfilter(I,PSF,'symmetric','conv');
V = 0*1E-2;
blurred_noisy = imnoise(blurred,'gaussian',0,V);

% Deconvolve
tic
Iout=deconvlucy(blurred_noisy,PSF,10);
toc

subplot(221)
imagesc(I);
axis equal tight
title('point sources');
colorbar

str=['$\sum \mathrm{{\bf I}_{i,j}}=' num2str(sum(sum(I))) '$'];
text(.02,.02,str,'units','normalized','color','w','verticalalignment','bottom',...
    'interpreter','latex','fontsize',12);


subplot(222);
imagesc(PSF);
colormap(purplemap);
axis equal tight
title('PSF (sum=1)');
colorbar


str=['$\sum \mathrm{{\bf PSF}_{i,j}}=' num2str(sum(sum(PSF))) ',~' ...
    '\sigma=' num2str(s) '$'];
text(.02,.02,str,'units','normalized','color','w','verticalalignment','bottom',...
    'interpreter','latex','fontsize',12);

subplot(223)
imagesc(blurred_noisy);
axis equal tight
title('blurred image w/ noise');
colorbar

str=['${\bf A}=\mathrm{{\bf PSF}}\star {\bf N}+\epsilon_{ij}$' newline ...
    '$\sum {\bf A}_{ij} = ' num2str(round(sum(sum(blurred_noisy)),1)) '$'];
text(.02,.02,str,'units','normalized','color','w','verticalalignment','bottom',...
    'interpreter','latex','fontsize',12);

subplot(224)
imagesc(Iout);
axis equal tight
title(['deconvolved']);
colorbar
% caxis([0 1]);

str=['${\widetilde {\bf I''}}=\widetilde{{\bf A}}/\widetilde{{\bf {\mathrm PSF}}}$' newline ...
    '$\sum{I''}=' num2str(round(sum(sum(Iout)),1)) '$'];
text(.02,.02,str,'units','normalized','color','w','verticalalignment','bottom',...
    'interpreter','latex','fontsize',12);

%% Deconvlution
mag=[82.6 83.2];         % Magnification x and y
pxsize=[16 16];          % Pixel size
pxsize_real=pxsize./mag; % um/pixel in real space
dPx=mean(pxsize_real);

% Get the sampling of the psf
dX=Xpsf(2)-Xpsf(1);
dY=Ypsf(2)-Ypsf(1);
dL=mean([dX dY]);
% I think this is units of nanometers

sc=2;
X=1:512;
Y=1:512;

% xL=[1 512];
% yL=[1 512];

Zsub=Z(yL(1):yL(2),xL(1):xL(2));
Xsub=xL(1):xL(2);
Ysub=yL(1):yL(2);

Zimg=imresize(Zsub,sc);
% Ximg=imresize(Xsub,sc);
% Yimg=imresize(Ysub,sc);



sfake=260*2; % gaussian radius of PSF fit

sfakePx=sfake*1E-3/(dPx/sc);

fakepsf = fspecial('gaussian',70,sfakePx);

hf4=figure;
clf
hf4.Color='w';
hf4.Name='Deconvolved Data';
colormap(purplemap);

subplot(221);
imagesc(Zimg);
axis equal tight
cc=[0 700];
colorbar
caxis(cc);
% xlim(xL);
% % ylim(yL);
drawnow;

subplot(222)
cla

imagesc(fakepsf);
axis equal tight

caxis([0 max(max(fakepsf))]);

colorbar
drawnow;

subplot(224)
tic
out = richlucy(Zimg,fakepsf,5);
toc
out2=imgaussfilt(out,0.5*sc);
imagesc(out2);
axis equal tight

caxis(cc);
colorbar
% xlim(xL);
% % ylim(yL);
disp('done');

%% Finding the lattice grid

