function [data] = processRawData(data,opts)

% This function takes the raw images and does processing on them.

if nargin == 1
    opts                    = struct;
    
    % Position Space
    opts.doSubtractBias     = 1;
    opts.doGaussFilter      = 0;
    opts.GaussFilterRadius  = 1;
    opts.doMask             = 0;
    opts.Mask               = ones(512,512);
    opts.doPSF              = 0;    
    opts.PSF                = [1.3 50 5];    
    opts.doRotate           = 0;
    opts.Theta              = 0;
        
    % Momentum Space
    opts.doFFT              = 1;
    opts.doMaskIR           = 1;
    opts.IRMaskRadius       = 0.01;
    opts.doFFTFilter        = 1;
    opts.FFTFilterRadius    = 1;    
end

if ~isfield(opts,'doRotate'); opts.doRotate = 1;end
    

for kk=1:length(data)

    fprintf([num2str(kk) ' of ' num2str(length(data)) ': '])
    t1=now;

    data(kk).Z = data(kk).RawImages;           
    data(kk).X = 1:size(data(kk).Z,2);
    data(kk).Y = 1:size(data(kk).Z,1); 
    
%% Remove extra exposure data
    if isfield(data(kk),'Flags') && isfield(data(kk).Flags,'qgm_doClearBufferExposure') 
        fprintf('removing extra exposure ...');
        if data(kk).Flags.qgm_doClearBufferEposure && size(data(kk).Z,3)>1
            data(kk).Z(:,:,1)=[];    
        end
    else
        if size(data(kk).Z,3)>1
            data(kk).Z(:,:,1)=[];   
        end
    end         
%% Subtract Digital Bias
    if opts.doSubtractBias  
        fprintf(' biasing ...');
        data(kk).Z = data(kk).Z - 200;
    end        
    
    %% Subtract Digital Bias
    if opts.doRotate  
        theta = opts.Theta;
        fprintf([' rotating (' num2str(theta) ' deg)...']);
        data(kk).Z = imrotate(data(kk).Z,theta,'bicubic');
    end        
%% Image Mask
    if opts.doMask
        fprintf('masking ...');
        for ii=1:size(data(kk).Z,3)  
            data(kk).Z(:,:,ii) = data(kk).Z(:,:,ii).*opts.Mask;        
        end
    end
%% Gaussian Fitler
    if opts.doGaussFilter  
        fprintf('smooth position ...');

        for ii=1:size(data(kk).Z,3)  
            data(kk).Z(:,:,ii) = imgaussfilt(data(kk).Z(:,:,ii),opts.GaussFilterRadius);
        end
    end    
    
%% Deconvolve PSF
    if opts.doPSF       
        fprintf('PSF ...');
        s       = opts.PSF(1);
        N       = opts.PSF(2);
        Niter   = opts.PSF(3);
        psf     = fspecial('gaussian',N,s);        
        for ii = 1:size(data(kk).Z,3)  
            data(kk).Z(:,:,ii) = deconvlucy(data(kk).Z(:,:,ii),psf,Niter);
        end       
    end

%% Compute FFT

    if opts.doFFT
        fprintf('FFT ...');
        % Compute FFT
        Nfft = 2^10;
        
        data(kk).Zf = zeros(Nfft,Nfft,size(data(kk).Z,3));
        data(kk).ZfNorm = zeros(Nfft,Nfft,size(data(kk).Z,3));
        data(kk).ZfPhase = zeros(Nfft,Nfft,size(data(kk).Z,3));

        for ii=1:size(data(kk).Z,3)
            zf = fft2(data(kk).Z(:,:,ii),Nfft,Nfft);
            zf = fftshift(zf);        
            data(kk).f = 1/2*linspace(-1,1,size(zf,1));       
            data(kk).Zf(:,:,ii) = zf;
            data(kk).ZfNorm(:,:,ii) = abs(zf);
            data(kk).ZfPhase(:,:,ii) = angle(zf);
        end   
    end
        
%% Mask IR
    if opts.doMaskIR && opts.doFFT
        fprintf('mask FFT ...');        
        kR = opts.IRMaskRadius;
        [fxx,fyy]=meshgrid(data(kk).f,data(kk).f);
        fmag = sqrt(fxx.^2+fyy.^2);
        iMask = fmag>kR;
        for ii=1:size(data(kk).Zf,3)
            data(kk).Zf(:,:,ii) = data(kk).Zf(:,:,ii).*iMask;
            data(kk).ZfNorm(:,:,ii) = data(kk).ZfNorm(:,:,ii).*iMask;
            data(kk).ZfPhase(:,:,ii) = data(kk).ZfPhase(:,:,ii).*iMask;
        end
    end
%% Filter FFT
    if opts.doFFT && opts.doFFTFilter
        fprintf('smooth FFT ...');
        for ii=1:size(data(kk).Zf,3)  
            data(kk).ZfNorm(:,:,ii) = imgaussfilt(data(kk).ZfNorm(:,:,ii),opts.FFTFilterRadius);
            data(kk).ZfPhase(:,:,ii) = imgaussfilt(data(kk).ZfPhase(:,:,ii),opts.FFTFilterRadius);
        end
    end
%% Finish it
    t2 = now;
    disp(['done (' num2str(round((t2-t1)*24*60*60,2)) ' s)']);     
    
        data(kk).X = 1:size(data(kk).Z,2);
    data(kk).Y = 1:size(data(kk).Z,1); 
end

end