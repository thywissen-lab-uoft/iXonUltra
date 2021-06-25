function ixondata = ixon_computeFFT(ixondata,opts)
disp('Computing 2d fft');

% Options for the 2D FFT
    if nargin==1
        opts=struct;
        opts.doSmooth=0;
        opts.smoothRadius=1;
        opts.fft_N=2^10;
    end
    
    % On NFFT= number of points in the FFT
    % It is typically useful to upscale the image to improve the smoothness
    % of the FFT (basically oversampling to get finer momentum resolution)
    
    % Sampling (1 pixel always)
    dL=1;  
    
    for kk=1:length(ixondata)
        fprintf([num2str(kk) ' of ' num2str(length(ixondata)) ' ... ']);
        % Compute the FFT
        fft_Z=fft2(ixondata(kk).Z,opts.fft_N,opts.fft_N);    
        
        % Renormalize the FFT
        fft_Z=fft_Z/(size(ixondata(kk).Z,1)^2);
        
        % Take absolute value
        fft_Z=abs(fft_Z);
        
        % Calculate the frequency vector
        fft_F=1/(2*dL)*linspace(-1,1,opts.fft_N);

        fft_Z=fftshift(fft_Z);  
        if opts.doSmooth
            fft_Z=imgaussfilt(fft_Z,opts.smoothRadius);
        end  
        
        % Assign output
        ixondata(kk).fft_Z=fft_Z;
        ixondata(kk).fft_F=fft_F;
        ixondata(kk).fft_N=opts.fft_N;
        disp('done');
    end

end

