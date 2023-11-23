function focus= ixon_focusStripe(data,stripe)
threshhold = 0.1;

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

z = imgaussfilt(z,1);

% z(z<10)=0;
[xx,yy]= meshgrid(x,y);


% Create phase and amplitude map
phiMap = stripe.PhaseMapFunc(L,theta,phi,xx,yy)+pi/2;
ampMap = stripe.EnvelopeFunc(A,xC,yC,s1,s2,theta,xx,yy)/A;

% Create a threshholding map to only look at statisfically significant
% regions
threshholdMap = ampMap>threshhold;


pL_N = floor(min(min(phiMap))/(2*pi))+1;
pH_N = floor(max(max(phiMap))/(2*pi))-1;
nVec = pL_N:pH_N;
stripe_boundary_lines = {};


x_coms = zeros(length(nVec),1);
y_coms = zeros(length(nVec),1);
z_stripes = zeros(length(y),length(x),length(nVec));
stripe_sum = zeros(length(nVec),1);
for n=1:length(nVec)
    nn=nVec(n);
    ii=[abs(phiMap-(nn*2*pi))<.3];
    stripe_boundary_lines{n}=polyfit(xx(ii),yy(ii),1);    
    i1 = (phiMap>=(nn*2*pi));
    i2 = (phiMap<=((nn+1)*2*pi));      

    stripe_map = i1.*i2;
    this_map=logical(stripe_map.*threshholdMap);

    % this_map=logical(stripe_map);

    z_stripes(:,:,n) = z.*this_map;
    stripe_sum(n) = sum(z_stripes(:,:,n),'all');
    x_coms(n) = sum(stripe_map.*xx,'all')/sum(stripe_map,'all');
    y_coms(n) = sum(stripe_map.*yy,'all')/sum(stripe_map,'all');
end

z0 = sum(z_stripes,'all');
%% Compute FFT of each stripe
Nfft = 2^10+1;
dX = x(2)-x(1);
f_max = 1/dX;
f = 1/2*linspace(-f_max,f_max,Nfft);    
[fxx,fyy]=meshgrid(f,f);
% Gaussian radius of PSF in real space in camera pixels
sPSF = 1.3161;

% Gaussian radius of ideal PSF in momentum space
qPSF = sqrt(1/(4*pi*sPSF^2))    
zf_psf = exp(-(fxx.^2+fyy.^2)/(2*qPSF^2));
zf_psf = sum(zf_psf.*zf_psf);

scores = zeros(length(nVec),1);
for n=1:size(z_stripes,3)
    figure(2001);

    if stripe_sum(n)<z0*0.05
        scores(n)=NaN;
        continue;
    end


    [sharpnessScore map]=MLVSharpnessMeasure(z_stripes(:,:,n));
    zf = fft2(z_stripes(:,:,n),Nfft,Nfft);
    zf = fftshift(zf); 
% zf = imgaussfilt(zf,2);
    zf = zf/stripe_sum(n);

    zfnorm=abs(zf);
    LPF = sqrt(fxx.^2+fyy.^2)<(2*qPSF);
    HPF = sqrt(fxx.^2+fyy.^2)>(.01);
    myfilt = LPF.*HPF;
    zf = zf.*myfilt;
    % zf = zf/stripe_sum(n);
    zfnorm2=abs(zf);
    
    % zfnorm2=zfnorm2/stripe_sum(n);
    scores(n) = norm(sum(zf.*zf_psf,'all'))
% scores(n)=sharpnessScore;
    subplot(size(z_stripes,3),2,2*(n-1)+1);
    imagesc(x,y,z_stripes(:,:,n));
        % axis equal tight

    subplot(size(z_stripes,3),2,2*n);
    imagesc(f,f,zfnorm2);
    caxis([0 1e-2])
    % % axis equal tigh
    
end




%%
focus = struct;
focus.scores = scores;
focus.x_coms = x_coms;
focus.y_coms = y_coms;
focus.sums = stripe_sum;

end

