function ixondata = ixon_computeFFT(ixondata,opts)

    if nargin==1
        opts=struct;
        opts.doSmooth=0;
        opts.smoothRadius=1;
    end

    for kk=1:length(ixondata)
        Zfft=fft2(ixondata(kk).Z);
        Zfft=abs(fftshift(Zfft));  
        if opts.doSmooth
            Zfft=imgaussfilt(Zfft,opts.smoothRadius);
        end        
        ixondata(kk).Zfft=Zfft;   
    end

end

