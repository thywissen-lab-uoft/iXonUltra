function [data] = ixonProcessImages(data,opts)

% This function takes the raw images and does processing on them.

if nargin == 1
    opts                    = struct;    
    % Position Space
    opts.doSubtractBias     = 1;
    opts.doSubtractBGD      = 0;
    opts.doScale            = 1;
    opts.ScaleFactor        = 2;
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

    data(kk).X = 1:size(data(kk).RawImages,2);
    data(kk).Y = 1:size(data(kk).RawImages,1);
    
    Z = data(kk).RawImages;
    
    size(Z)
    
%% Subtract Digital Bias
    if opts.doSubtractBias  
        fprintf(' biasing ...');
        Z = Z - 200;
    end      
    
%% Remove wipe pics

if isfield(data.Params,'qgm_MultiExposures') && isfield(data.Params,'qgm_MultiPiezos')
    wipePics = isnan([data.Params.qgm_MultiExposures]);
    Z(:,:,wipePics)=[];
    L = size(Z,3);
    if isfield(data(kk).Flags,'lattice_fluor_bkgd') && mod(L,2)==0
        Zme = zeros(size(Z,1),size(Z,2),L/2);
        for n=1:L/2
            Zme(:,:,n) = Z(:,:,n) -  Z(:,:,n+L/2);
        end
        Z=Zme;
    else
        Z = Zraw;
    end
    
else
    if isfield(data(kk),'Flags') && isfield(data(kk).Flags,'lattice_ClearCCD_IxonTrigger') 
        fprintf('removing extra exposure ...');
        if data(kk).Flags.lattice_ClearCCD_IxonTrigger && size(data(kk).Z,3)>1
            data(kk).Z(:,:,1)=[];    
        end
    else
         if size(data(kk).Z,3)>1
             data(kk).Z(:,:,1)=[];   
         end
    end   
end

data(kk).Z = Z;
    %%
    
%     data(kk).Z = data(kk).RawImages;           



%% Raw Data
    data(kk).Zraw = data(kk).Z;    
    
%% Remove background exposure data
    if isfield(data(kk),'Flags') && isfield(data(kk).Flags,'lattice_fluor_bkgd') 
        if data(kk).Flags.lattice_fluor_bkgd
            data(kk).BG = data(kk).Z(:,:,end);
            if opts.doSubtractBG
                fprintf('removing background exposure ...');
                numImag = size(data(kk).Z,3);
                for ii=1:(numImag-1)
                    data(kk).Z(:,:,ii)=data(kk).Z(:,:,ii)-data(kk).BG;
                end
                data(kk).Z(:,:,end) = [];
            end
        end
    end 

%% Image Mask
    if opts.doMask
        fprintf('masking ...');
        for ii=1:size(data(kk).Z,3)  
            data(kk).Z(:,:,ii) = data(kk).Z(:,:,ii).*opts.Mask;        
        end
    end

%% Scale Image
    if isfield(opts,'doScale') && opts.doScale         
       data(kk).X = linspace(1,length(data(kk).X),length(data(kk).X)*opts.ScaleFactor);
       data(kk).Y = linspace(1,length(data(kk).Y),length(data(kk).Y)*opts.ScaleFactor);
       data(kk).Z = imresize(data(kk).Z,opts.ScaleFactor)/(opts.ScaleFactor)^2;
    end    
%% Rotate Image
    if opts.doRotate  
        theta = opts.Theta;
        fprintf([' rotating (' num2str(theta) ' deg)...']);
        data(kk).Z = imrotate(data(kk).Z,theta,'bicubic','crop');
        data(kk).RotationMask = imrotate(ones(size(data(kk).Z,1),size(data(kk).Z,1)),theta,'bicubic','crop');
        data(kk).RotationMask = logical(round(data(kk).RotationMask));
    end           
%% No Filter
    data(kk).ZNoFilter = data(kk).Z; 

%% Gaussian Fitler
    if opts.doGaussFilter  
        fprintf('smooth position ...');
        for ii=1:size(data(kk).Z,3)  
            data(kk).Z(:,:,ii) = imgaussfilt(data(kk).Z(:,:,ii),opts.GaussFilterRadius);
        end
    end        
%% Deconvolve Point Spread Function
    if opts.doPSF       
        fprintf('PSF ...');
        s       = opts.PSF(1);
        if opts.doScale
           s = s*opts.ScaleFactor; 
        end
        N       = opts.PSF(2);
        Niter   = opts.PSF(3);
        psf     = fspecial('gaussian',N,s);           
        for ii = 1:size(data(kk).Z,3)  
            data(kk).Z(:,:,ii) = deconvlucy(data(kk).Z(:,:,ii),psf,Niter); 
        end        
    end  
%% Fast Fourier Transform (FFT)
    if opts.doFFT
        fprintf('FFT ...');
        % Compute FFT
        Nfft = 2^10+1;        
        data(kk).Zf = zeros(Nfft,Nfft,size(data(kk).Z,3));
        data(kk).ZfNorm = zeros(Nfft,Nfft,size(data(kk).Z,3));
        data(kk).ZfPhase = zeros(Nfft,Nfft,size(data(kk).Z,3));
        for ii=1:size(data(kk).Z,3)
            zf = fft2(data(kk).Z(:,:,ii),Nfft,Nfft);
            zf = fftshift(zf);                    
            dX = data(kk).X(2)-data(kk).X(1);
            f_max = 1/dX;
            data(kk).f = 1/2*linspace(-f_max,f_max,size(zf,1));              
            data(kk).Zf(:,:,ii) = zf;
            data(kk).ZfNorm(:,:,ii) = abs(zf);
            data(kk).ZfPhase(:,:,ii) = angle(zf);
        end   
    end        
%% Mask IR in FFT
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
    data(kk).ProcessOptions = opts;
%         data(kk).X = 1:size(data(kk).Z,2);
%     data(kk).Y = 1:size(data(kk).Z,1); 
end

end