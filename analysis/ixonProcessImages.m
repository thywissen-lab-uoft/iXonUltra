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
    if ~isfield(data(kk),'Flags')    
        data(kk).Flags = struct;
    end


    tic
    fprintf([num2str(kk) ' of ' num2str(length(data)) ': '])
    data(kk).X = 1:size(data(kk).RawImages,2);
    data(kk).Y = 1:size(data(kk).RawImages,1);   
    
    Z = data(kk).RawImages;   
    
%% Subtract Digital Bias
% Subtract digital bias of 200 counts
    if opts.doSubtractBias  
        fprintf(' biasing ...');
        Z = Z - 200;
    end      
    
%% Remove wipe pics

    if isfield(data(kk).Params,'qgm_MultiExposures') && isfield(data(kk).Params,'qgm_MultiPiezos')
        wipePics = isnan([data(kk).Params.qgm_MultiExposures]);
        Z(:,:,wipePics)=[];      
    elseif isfield(data(kk).Flags,'lattice_ClearCCD_IxonTrigger')  && ...
            data(kk).Flags.lattice_ClearCCD_IxonTrigger
            if data(kk).Flags.lattice_ClearCCD_IxonTrigger && size(Z,3)>1
                Z(:,:,1)=[];    
            end
    elseif ~ isfield(data(kk).Flags,'lattice_ClearCCD_IxonTrigger') && ...
        size(data(kk).Z,3)>1        
            Z(:,:,1)=[]; 
    end

    data(kk).Z = Z;

%% Raw Data
    data(kk).Zraw = data(kk).Z;    
    
%% Remove background exposure data

if opts.doSubtractBG 
    if isfield(data(kk),'lattice_fluor_bkgd') && data(kk).lattice_fluor_bkgd == 0
        warning('Trying to remove bkgd image but flag appears to be deactivated.  If the flag is depreciated, please remove this check.');
    end
    Z = data(kk).Z;
    L = size(Z,3);
    if mod(L,2)~=0
        error('Cannot subtract background with odd number of images');
    end
    Z_bg_removed = zeros(size(Z,1),size(Z,2),L/2);
    for n=1:L/2
        Z_bg_removed(:,:,n) = Z(:,:,n) -  Z(:,:,n+L/2);
    end
    Z=Z_bg_removed;
    data(kk).Z = Z_bg_removed;

    % Measure the noise of the subtraction
    for n = 1:L/2
        Zthis = Z(:,:,n);
        Zneg = Zthis(Zthis<0);
        [N,edges] = histcounts(Zneg);
        % iZero = find(edges>=0,1);
        centers = (edges(1:end-1) + edges(2:end))/2;

        vals = cumsum(N,'reverse')/sum(N);
        ind = find(vals<=.5,1);
        thresh = -centers(ind);
        data(kk).NoiseEstimation(n) = thresh;      
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
        data(kk).X = linspace(1,data(kk).X(end),length(data(kk).X)*opts.ScaleFactor);
       data(kk).Y = linspace(1,data(kk).Y(end),length(data(kk).Y)*opts.ScaleFactor);
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
            data(kk).Z(:,:,ii) = imgaussfilt(data(kk).Z(:,:,ii),...
                opts.GaussFilterRadius,'FilterDomain','frequency');
        end
    end        
%% Deconvolve Point Spread Function
    if opts.doPSF       
        fprintf('PSF ...');
        s       = opts.PSF(1);    
        if opts.doScale;s = s*opts.ScaleFactor;end
        N       = opts.PSF(2);
        Niter   = opts.PSF(3);
        psf     = fspecial('gaussian',N,s);  
        for ii = 1:size(data(kk).Z,3) 
            if opts.doSubtractBG 
                Nnoise=data(kk).NoiseEstimation(ii);
            else
                Nnoise = 50;
            end           
            noise_variance = Nnoise^2;

            % data(kk).Z(:,:,ii)=data(kk).Z(:,:,ii)+Nnoise;
            % noise_variance=noise_variance*2;

            data(kk).Z(:,:,ii) = deconvlucy(data(kk).Z(:,:,ii),...
                psf,Niter,0,1,noise_variance);      
            % data(kk).Z(:,:,ii)=data(kk).Z(:,:,ii)-Nnoise;

        end        
    end  
%% Fast Fourier Transform (FFT)
    if opts.doFFT
        fprintf('FFT ...');
        % Compute FFT
        Nfft = 2^10+1;        

        src = data(kk).Z;
        l = size(src,3);

        data(kk).Zf = zeros(Nfft,Nfft,l);
        data(kk).ZfNorm = zeros(Nfft,Nfft,l);
        data(kk).ZfPhase = zeros(Nfft,Nfft,l);
        for ii=1:l
            zf = fft2(src(:,:,ii),Nfft,Nfft);
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
    t2 = toc;
    disp(['done (' num2str(round(t2,2)) ' s)']);         
    data(kk).ProcessOptions = opts;
%         data(kk).X = 1:size(data(kk).Z,2);
%     data(kk).Y = 1:size(data(kk).Z,1); 
end

end