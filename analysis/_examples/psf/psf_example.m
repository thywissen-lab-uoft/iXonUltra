function psf_example
% PSF and Deconvoution
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

%% Creat the PSF and fake data
% Create a PSF
N=70; % size of PSF (NxN)
s=5;  % gaussian radius (pixels)
PSF = fspecial('gaussian',N,s);
nLucy=5;

PSF=PSF;

%% Prepare Image

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


%% Deconvolve
tic
Iout=deconvlucy(blurred_noisy,PSF,nLucy);
toc
%%

hf5=figure;
clf
hf5.Color='w';
hf5.Name='PSF Example';
hf5.Position=[100 100 1000 800];

% Plot the source data
subplot(221)
imagesc(I);
axis equal tight
title('point sources');
colorbar

str=['$\sum \mathrm{{\bf I}_{i,j}}=' num2str(sum(sum(I))) '$'];
text(.02,.02,str,'units','normalized','color','w','verticalalignment','bottom',...
    'interpreter','latex','fontsize',12);

% Plot the PSF
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

% Plot the blurred image
subplot(223)
imagesc(blurred_noisy);
axis equal tight
title('blurred image w/ noise');
colorbar

str=['${\bf A}=\mathrm{{\bf PSF}}\star {\bf N}+\epsilon_{ij}$' newline ...
    '$\sum {\bf A}_{ij} = ' num2str(round(sum(sum(blurred_noisy)),1)) '$'];
text(.02,.02,str,'units','normalized','color','w','verticalalignment','bottom',...
    'interpreter','latex','fontsize',12);

% Plot the deconvolved data
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


end

