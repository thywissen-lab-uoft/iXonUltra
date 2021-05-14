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
xL=[310 400];
yL=[340 410];
cL1=[0 200];
cL2=[0 100];

hf=figure(1);
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
%    xlim(xL);
%    ylim(yL);
end
colormap(purplemap);



%%
% Load point spread function
psf='Y:\Lattice Experiment Programs\Lattice Image Analysis\GUI 1-2-2\Jan_29_PSF.mat';
psf=load(psf);
psf=psf.psf;

hf2=figure(2);
hf2.Name='Measured Point Spread Function';
clf
Zpsf=psf.psf;
Zpsf=Zpsf/max(max(Zpsf));
Xpsf=psf.x(1,:);
Ypsf=psf.y(:,1);
imagesc(Xpsf,Ypsf,Zpsf)

axis equal tight

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

hf5=figure(5);
clf
hf5.Color='w';
hf5.Name='PSF Example';

N=1000; % size of PSF (NxN)
s=100;  % gaussian radius (pixels)
PSF = fspecial('gaussian',1000,100);

supblot(221)
blurred = imfilter(I,PSF,'symmetric','conv');

subplot(222);
imagesc(PSF);
colormap(purplemap);



%% Deconvlution


hf4=figure(4);
clf
hf4.Color='w';
hf4.Name='Deconvolved Data';
out = richlucy(Z,1,5);

imagesc(out)

colormap(purplemap);
caxis(cL1);


