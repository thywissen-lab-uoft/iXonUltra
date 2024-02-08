function ixondata=sharpenPSF(ixondata)
disp('Sharpening with PSF');

doDebug=1;
%% Image Parameters
% These will be the assumed image parameters

mag=[82.6 83.2];         % Magnification x and y
pxsize=[16 16];          % Pixel size
pxsize_real=pxsize./mag; % um/pixel in real space
dPx=1e3*mean(pxsize_real);   % Size of pixel in nm

%% Lucy Number 
nLucy=4;

%% Rescale Settings
% It is useful to rescale the image so that the PSF is super resolved.
% However, this increas the computation time
sc=3;
dPx=dPx/sc;



%% Gaussian PSF
% Create the PSF and rescale it so that each pixel is the same as the
% scaled data image.

% Gaussian radius of PSF in nm (Graham's Thesis);
sPSF=260;
% sPSF=200;

% Choose the size of the PSF be 6*sigma
Npx=sPSF*6/dPx;
Npx=round(Npx);

% Make it an odd matrix so that zero is included.
Npx=(1-mod(Npx,2))+Npx;

gauss_psf = fspecial('gaussian',Npx,sPSF/dPx);

xpsf=dPx*(1:Npx);
xpsf=xpsf-mean(xpsf); 

%% Show PSf

hF=figure;
hF.Color='w';

imagesc(xpsf,xpsf,gauss_psf);
xlabel('nm');
ylabel('nm');
colormap(purplemap);

title('PSF');

%% Sharpen with PSF

psf=gauss_psf;
if doDebug
    hF_debug=figure;
    hF_debug.Position=[100 100 1200 800];
    hF_debug.Color='w';
end

for kk=1:length(ixondata)
    fprintf([num2str(kk) ' of ' num2str(length(ixondata)) ' ... ']);      
    Z=ixondata(kk).Z;
    X=1:size(Z,2);X=X';
    Y=1:size(Z,1);Y=Y';
    
    if doDebug
       figure(hF_debug);
       colormap(purplemap);
       subplot(121);
       imagesc(X,Y,Z);
      caxis([0 1500]);
             axis equal tight

    end    
    
    % Resize the data
    Z=imresize(Z,sc);
    X=imresize(X,[length(X)*sc 1]);
    Y=imresize(Y,[length(Y)*sc 1]);
    
    % Perform the deconvolution
    Zout = richlucy(Z,psf,nLucy);
    
    % Resize the data back down
    Zout = imresize(Zout,1/sc);
    
    Zout = imgaussfilt(Zout,.75);
    
    if doDebug
       figure(hF_debug);
       subplot(122);
       imagesc(X,Y,Zout);
      caxis([0 1500]);
      axis equal tight
    end
    
    
    ixondata(kk).Z=Zout;
    ixondata(kk).Zold=Z;
    disp('done');
    
    if doDebug
        waitforbuttonpress
    end
end
disp('hi');

end

