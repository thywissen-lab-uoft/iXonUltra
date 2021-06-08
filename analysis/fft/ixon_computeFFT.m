function ixondata = ixon_computeFFT(ixondata,opts)

    if nargin==1
        opts=struct;
    end

    for kk=1:length(ixondata)
        Zfft=fft2(Z);
        Zfft=abs(fftshfit(Y));
        ixondata(kk).Zfft=Zfft;    
    end

end

