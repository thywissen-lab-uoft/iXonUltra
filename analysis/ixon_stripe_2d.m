%% 2D Stripe Analysis
% This analzyes the stripes for take from the ixon camera.  This is a 2D
% fit over the entire cloud.  This fits a 2D gaussian modulated by a sine
% wave at a particular angle.  This is pariticularly useful to fit the
% angular dependence of the data. 
%
% The input fit parameters are specified in the options structure.


stripe_2d_opts=struct;

stripe_2d_opts.xUnit=ixon_unit;

stripe_2d_opts.ShimFit=0;
stripe_2d_opts.Theta=[10 100]; % Specify the domain (MUST BE 180 DEGREES)

stripe_2d_opts.saveAnimation=1;        % save the animation?
stripe_2d_opts.StartDelay=1;
stripe_2d_opts.MidDelay=.25;
stripe_2d_opts.EndDelay=1;
stripe_2d_opts.CLim = [0 100];

stripe_2d_opts.ROI = NaN;
% stripe_2d_opts.ROI = [100 450 25 425];

stripe_2d_opts.Threshhhold = 20; % Ignore pixels below this threshhhold
stripe_2d_opts.Threshhhold = NaN;

stripe_2d_opts.ConstrainWavelength = NaN;
% stripe_2d_opts.ConstrainWavelength = 72;

stripe_2d_opts.ConstrainAngle = NaN;
%  stripe_2d_opts.ConstrainAngle = 88;

ixondata=flip(ixondata);
if ixon_doAnalyzeStripes2D
    [hF_stripe_2d,stripe_data2d]=analyzeStripes2(ixondata,ixon_xVar,stripe_2d_opts);

    if ixon_doSave
        ixon_saveFigure(ixondata,hF_stripe_2d,'ixon_field_stripe_2d');        
    end
    
    outdata.stripe_data2d=stripe_data2d;   
end

%% Stripe Analysis FFT
ixon_doAnalyzeStripes2D_Focusing=0;
if ixon_doAnalyzeStripes2D_Focusing
    for kk=90:90
        Z = ixondata(kk).Z;
        x=1:size(Z,2);x=x';
        y=1:size(Z,1);y=y';    
        [xx,yy]=meshgrid(x,y);
        theta = stripe_data2d.Theta(kk,1);
        phi = stripe_data2d.Phi(kk,1);
        L = stripe_data2d.Wavelength(kk,1);
        
        phaseMap = 2*pi/L*(cosd(theta)*xx+sind(theta)*yy)+phi+pi/2;
        
        pL_N = floor(min(min(phaseMap))/(2*pi))+1;
        pH_N = floor(max(max(phaseMap))/(2*pi))-1;
        sVec = [];
        Nvec = [];
        nVec = pL_N:pH_N;
        for nn=pL_N:pH_N
            i1 = (phaseMap>=(nn*2*pi));
            i2 = (phaseMap<=((nn+1)*2*pi));            
            this_map=logical(i1.*i2);
            
            Z_onestripe = Z.*this_map;
                    Nfft = 2^10;



            zf = fft2(Z_onestripe,Nfft,Nfft);
            zf = fftshift(zf);       
            
            f = 1/2*linspace(-1,1,size(zf,1));    
            zfnorm = abs(zf);
            
            zfnorm=zfnorm/sum(sum(zfnorm));
            
            
            
            
            [f1,f2] = meshgrid(f,f);
            
%             fmag = sqrt(f1.^2+f2.^2);
%             ifilt = [fmag>.1];
%             
%             zfnorm = ifilt.*zfnorm;
%                         zfnorm=zfnorm/sum(sum(zfnorm));

            
            v1 = sum(sum(zfnorm.*f1.^2));
            v2 = sum(sum(zfnorm.*f1.^2));            
            s = sqrt(v1+v2);            
            sVec(end+1)=s;
            Nvec(end+1)=sum(sum(Z_onestripe));
        end
        
        hF = figure(917);
        clf
        hF.Color='w';
        
        subplot(221);
        imagesc(x,y,Z);
        set(gca,'ydir','normal');
        caxis([0 100])
        
        subplot(222);
        [~,ind]=min(sVec);
        nn=nVec(ind);
                 i1 = (phaseMap>=(nn*2*pi));
            i2 = (phaseMap<=((nn+1)*2*pi));            
            this_map=logical(i1.*i2);
            
            Z_onestripe = Z.*this_map;
        
        imagesc(x,y,Z_onestripe);
        set(gca,'ydir','normal');
                caxis([0 100])

        subplot(223);
        yyaxis left
        plot(nVec,sVec,'o')
        ylabel('fft second moment');

        yyaxis right
        plot(nVec,Nvec,'o')
        ylabel('counts');
        xlabel('stripe index');
   
    end
end