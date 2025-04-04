function [data] = ixon_ProcessImages(data,opts)

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
    %     opts.Mask               = ones(512,512);
    load('C:\Users\Sephora\Documents\GitHub\iXonUltra\analysis\ixon_mask.mat');
    opts.Mask = BW;
    opts.doPSF              = 0;    
    opts.PSF                = [1.3 50 5];    
    opts.DetectNoise        = 1;

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
        Z = Z - 200; % This is a digital biasing (ie. from computer)
    end      
    
%% Remove wipe pics

if size(Z,3)>1
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
    Z_bg = zeros(size(Z,1),size(Z,2),L/2);

    R = 256 + [-100:100];


    for n=1:L/2
        Z_bg = Z(:,:,n+L/2);
        Z_bg_blur = imgaussfilt(Z_bg,2);
        dZ_bg = Z_bg-Z_bg_blur;
        pd = fitdist(dZ_bg(:),'normal');
        data(kk).NoiseEstimation(n) = pd.sigma;
        Z_bg_removed(:,:,n) = Z(:,:,n) - Z_bg_blur;
    end
    Z=Z_bg_removed;
    data(kk).Z = Z_bg_removed;



    % % Measure the noise of the subtraction
    % for n = 1:L/2
    %     try
    % 
    % 
    % 
    %         Zthis = Z(:,:,n);
    % 
    %         d=load('ixon_mask.mat');
    %         BW = d.BW;
    %         BW=double(BW);
    %         BW(BW==0)=nan;
    %         Zb = Zthis.*BW;
    %         Zc = Zb(:);
    %         Zc(isnan(Zc)) = [];            
    %         Zd = Zc;
    %         Zd(Zd>0)=[];
    %         Zd = [Zd; -Zd];            
    %         s = std(Zd);            
    %         data(kk).NoiseEstimation(n) = s ;
    % 
    %         Z_bg
    %     catch ME
    %         data(kk).NoiseEstimation(n) = 0;
    %     end
    % end
else
    data(kk).NoiseEstimation = 0;
end
  
%% Image Mask
    if opts.doMask
        fprintf('masking ...');
        for ii=1:size(data(kk).Z,3)  
            data(kk).Z(:,:,ii) = data(kk).Z(:,:,ii).*opts.Mask;        
        end
    end
    
%% No Filter
    data(kk).ZNoFilter = data(kk).Z; 

    %% Deconvolve Point Spread Function
    if opts.doPSF       
        fprintf('PSF ...');
        s       = opts.PSF(1);    
%         if opts.doScale;s = s*opts.ScaleFactor;end
        N       = opts.PSF(2);
        Niter   = opts.PSF(3);
        psf     = fspecial('gaussian',N,s);  
        for ii = 1:size(data(kk).Z,3) 
            
            % if  opts.DetectNoise  && opts.doSubtractBG
            %     Nnoise=data(kk).NoiseEstimation(ii);
            % 
            % else
            %     Nnoise = opts.Noise;
            %     data(kk).NoiseEstimation(ii) = Nnoise;
            % end            
            % noise_variance = Nnoise^2;
            % Zpre = data(kk).Z(:,:,ii);
            % Zpre(Zpre<=0)=0;
            % Zpre(Zpre<=(Nnoise))=0;          
            % 
            % Zsharp = deconvlucy(Zpre,...
            %     psf,Niter,0,1,noise_variance);             
            % data(kk).Z(:,:,ii) =    Zsharp;

            Zpre = data(kk).Z(:,:,ii);
            Nsigma_offset = 2;
            NoiseSigma = data(kk).NoiseEstimation(n);
            NoiseVariance = NoiseSigma^2;
            offset = Nsigma_offset*NoiseSigma;
            Zsharp = deconvlucy(Zpre+offset,...
                psf,Niter,0,1,NoiseVariance+offset)-offset;      
            data(kk).Z(:,:,ii) =    Zsharp;
        end        
    end      
  

%% Scale Image
    if isfield(opts,'doScale') && opts.doScale         
    
        data(kk).X = linspace(1,data(kk).X(end),length(data(kk).X)*opts.ScaleFactor);
       data(kk).Y = linspace(1,data(kk).Y(end),length(data(kk).Y)*opts.ScaleFactor);
       data(kk).Z = imresize(data(kk).Z,opts.ScaleFactor)/(opts.ScaleFactor)^2;
       
       data(kk).ZNoFilter = imresize(data(kk).ZNoFilter,opts.ScaleFactor)/(opts.ScaleFactor)^2;

    end    
%% Rotate Image
    if opts.doRotate  
        theta = opts.Theta;
        fprintf([' rotating (' num2str(theta) ' deg)...']);
        data(kk).Z = imrotate(data(kk).Z,theta,'bicubic','crop');
        data(kk).ZNoFilter = imrotate(data(kk).ZNoFilter,theta,'bicubic','crop');

        data(kk).RotationMask = imrotate(ones(size(data(kk).Z,1),size(data(kk).Z,2)),theta,'bicubic','crop');
        data(kk).RotationMask = logical(round(data(kk).RotationMask));        
    end           

%% Gaussian Fitler
    if opts.doGaussFilter  
        fprintf('smooth position ...');
        for ii=1:size(data(kk).Z,3)  
            data(kk).Z(:,:,ii) = imgaussfilt(data(kk).Z(:,:,ii),...
                opts.GaussFilterRadius,'FilterDomain','frequency');
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