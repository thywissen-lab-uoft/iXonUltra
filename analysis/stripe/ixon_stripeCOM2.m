function [ixondata,stripes,qgmdata_stripe] = ixon_stripeCOM2(ixondata,stripes,opts)
threshold_stripeROI = 0.1;

for jj=1:length(ixondata)


    %% Find stripe ROIs %%
    data = ixondata(jj);
    stripe = stripes(jj);

    % Fit Parameters
    theta = stripe.theta;
    phi = stripe.phi;
    L = stripe.L;
    xC = stripe.xC;
    yC = stripe.yC;
    s1 = stripe.s1;
    s2 = stripe.s2;
    A = stripe.A;
    
    % Grab Data
    x = data.X;
    y = data.Y;
    z = data.ZNoFilter;
    % z = data.Z;
    
    % z(z<=10)=0;
    % z=imgaussfilt(z,1);
    
    [xx,yy]= meshgrid(x,y);
    
    % Fourier Filtering
    Nfft = 2^10;
    sPSF = 1.3161;% Gaussian radius of PSF in real space in camera pixels
    qPSF = sqrt(1/(4*pi*sPSF^2));% Gaussian radius of ideal PSF in momentum space
    dX = x(2)-x(1);
    f_max = 1/dX;
    f = 1/2*linspace(-f_max,f_max,Nfft);    
    [fxx,fyy]=meshgrid(f,f);
    zf = fftshift(fft2(z,Nfft,Nfft));
    zf_psf = exp(-(fxx.^2+fyy.^2)/(2*qPSF.^2));
    
    
    lpf = ones(length(fyy),length(fyy));
    hpf = ones(length(fyy),length(fyy));
    
    % lpf = exp(-(fxx.^2+fyy.^2)/(2*qPSF.^2));
    % hpf = 1-exp(-(fxx.^2+fyy.^2)/(2*.01));
    
    
    zf = zf.*hpf;
    zf = zf.*lpf;
    
    
    z = abs(ifft2(zf,length(y),length(x)));
         % z=imgaussfilt(z,1);   
        
    % Create phase and amplitude map
    phiMap = stripe.PhaseMapFunc(L,theta,phi,xx,yy)+pi/2;
    ampMap = stripe.EnvelopeFunc(A,xC,yC,s1,s2,theta,xx,yy)/A;
    
    % Create a threshholding map to only look at statisfically significant
    % regions
    threshholdMap = ampMap>threshold_stripeROI;
    
    
    pL_N = floor(min(min(phiMap))/(2*pi))+1;
    pH_N = floor(max(max(phiMap))/(2*pi))-1;
    nVec = pL_N:pH_N;

    stripe_boundary_lines = {};
    
    
    x_coms = [];
    y_coms = [];
    z_stripes = zeros(length(y),length(x),length(nVec));
    stripe_sum = [];
    
    xBs=zeros(length(nVec),2);
    yBs = zeros(length(nVec),2);
    
    z_stripes={};
    x_stripes ={};
    y_stripes = {};

%     ROI = zeros(1,4,length(nVec));

    ROI = [];

    for n=1:length(nVec)
        nn=nVec(n);
        ii=[abs(phiMap-(nn*2*pi))<.3];
        stripe_boundary_lines{n}=polyfit(xx(ii),yy(ii),1);    
        i1 = (phiMap>=(nn*2*pi));
        i2 = (phiMap<=((nn+1)*2*pi));      
    
        % Create mask for this stripe
        stripe_map = i1.*i2;
        this_map=logical(stripe_map.*threshholdMap);    
        
        if sum(this_map,'all')==0

%             ROI(:,:,n) = [];
            continue;
        end
        
        indR = 1:size(this_map,1);
        indC = 1:size(this_map,2);
        
        [cc,rr] = meshgrid(indC,indR);
        
        cAll = cc(this_map); 
        rAll = rr(this_map);
        
        r1 = min(rAll(:));
        r2 = max(rAll(:));
        
        c1 = min(cAll(:));
        c2 = max(cAll(:));
        
        
        zThis = z(r1:r2,c1:c2);
        
        if size(zThis,1)<50 || size(zThis,2)<50 

           continue 
        end
        
        xThis = x(c1:c2);
        yThis = y(r1:r2);

        PSF = fspecial('gaussian',20,5);
        % zThis = edgetaper(zThis,PSF);
    
        
        z_stripes{end+1} = zThis;
        x_stripes{end+1} =xThis;
        y_stripes{end+1} = yThis;
        
        [xx2,yy2]=meshgrid(xThis,yThis);
        
        stripe_sum(end+1) = sum(zThis,'all');
    
        x_coms(end+1) = sum(xx2.*zThis,'all')/sum(zThis,'all');
        y_coms(end+1) = sum(yy2.*zThis,'all')/sum(zThis,'all');

        if isempty(ROI)
            ROI(1,:,1) = round([xThis(1) xThis(end) yThis(1) yThis(end)]);
        else
            ROI(1,:,end+1) = round([xThis(1) xThis(end) yThis(1) yThis(end)]);
        end
        
    
    end

%     [y_sort,sort_i] = sort(abs(y_coms-yC));
% 
%     x_coms = x_coms(sort_i);
%     y_coms = y_coms(sort_i);
%     ROI = ROI(:,:,sort_i);
%     stripes(jj).ind  = sort_i;

    stripes(jj).ROI = ROI;
    stripes(jj).xCOM = x_coms;
    stripes(jj).yCOM = y_coms;

   
    
    z0 = sum(stripe_sum,'all');


        %% Bin into Lattice Sites
% Having calculated the lattice spacing and phase. Bin all counts into a
% lattice site.


    nStripes = length(stripes(jj).xCOM);


    for ii=1:nStripes
        for kk=1:1
            opts = struct;
            a1 = ixondata(jj).LatticePhase(kk).a1;
            a2 = ixondata(jj).LatticePhase(kk).a2;                        
            p1 = ixondata(jj).LatticePhase(kk).p1;
            p2 = ixondata(jj).LatticePhase(kk).p2;        
            opts.ScaleFactor = 3;    
            opts.a1 = a1;
            opts.a2 = a2;
            opts.p1 = p1;
            opts.p2 = p2;  

            if isfield(ixondata(jj),'RotationMask')
               opts.Mask =  ixondata(jj).RotationMask;
            end

            ROI=stripes(jj).ROI(:,:,ii);

            ix_1 = find(ixondata(jj).X>=ROI(1),1);
            ix_2 = find(ixondata(jj).X>=ROI(2),1);
            iy_1 = find(ixondata(jj).Y>=ROI(3),1);
            iy_2 = find(ixondata(n).Y>=ROI(4),1);
            x = ixondata(jj).X(ix_1:ix_2);
            y = ixondata(jj).Y(iy_1:iy_2);   
            z = ixondata(jj).Z(iy_1:iy_2,ix_1:ix_2,kk);    
            tic;
            fprintf(['(' num2str(ii) '/' num2str(nStripes) ...
                ') binning stripe into lattice ...']);   
            qgmdata_stripe(jj).LatticeBin(ii) = binLattice(x,y,z,opts); 
            t2=toc;
            disp(['done (' num2str(t2,3) ' sec.)']);   

       

        end    
    end
    

%%  Digitize into lattice sites


    threshold_dig = 1500;


    nStripes = length(stripes(jj).xCOM);

    track_good = zeros(1,nStripes);

    for ii = 1:nStripes

        LatticeDig = struct;
  
        x = qgmdata_stripe(jj).LatticeBin(ii).n1;
        y = qgmdata_stripe(jj).LatticeBin(ii).n2;
        Zdig = qgmdata_stripe(jj).LatticeBin(ii).Zbin>=threshold_dig; 
        Natoms = sum(sum(Zdig));        % Total number of atoms




        zY=sum(Zdig,2)';zY = zY/sum(zY);
        zX=sum(Zdig,1); zX = zX/sum(zX);

        % Calculate center of mass
        Xc=sum(zX.*x);
        Yc=sum(zY.*y);          

        % Calculate central second moment/variance and the standard
        % deviation
        X2=sum(zX.*(x-Xc).^2); % x variance
        Xs=sqrt(X2); % standard deviation X
        Y2=sum(zY.*(y-Yc).^2); % x variance
        Ys=sqrt(Y2); % standard deviation Y               

        LatticeDig.Zdig = Zdig;
        LatticeDig.Natoms = Natoms;
        LatticeDig.n1 = x;
        LatticeDig.n2 = y;
        LatticeDig.Xc = Xc;
        LatticeDig.Yc = Yc;
        LatticeDig.Xs = Xs;
        LatticeDig.Ys = Ys;  


        
        if Natoms>2

            track_good(ii) = 1;
            qgmdata_stripe(jj).LatticeDig(sum(track_good)) = LatticeDig;

        end

            

    end

    fprintf('%.0f good stripes in run %.0f  \n',[sum(track_good) jj]);
    
    %clear info of bad stripes
    stripes(jj).ROI = stripes(jj).ROI(:,:,find(track_good));
    stripes(jj).xCOM = stripes(jj).xCOM(find(track_good));
    stripes(jj).yCOM = stripes(jj).yCOM(find(track_good));

    x_stripes = x_stripes(find(track_good));
    y_stripes = y_stripes(find(track_good));
    z_stripes = z_stripes(find(track_good));

        %% Show all Stripes
    
    hFme = figure(2001);
    clf
    
    for kk=1:length(stripes(jj).xCOM)

       ld = qgmdata_stripe(jj).LatticeDig(kk);

       subplot(length(z_stripes),2,2*kk-1);
       imagesc(x_stripes{kk},y_stripes{kk},z_stripes{kk}); 
       axis equal tight
       title('Stripe '+string(kk)); 

       subplot(length(z_stripes),2,2*kk);
       imagesc(ld.n1,ld.n2,ld.Zdig); 
       axis equal tight
       title('Stripe '+string(kk));


      
    end

    keyboard;

end   

end
