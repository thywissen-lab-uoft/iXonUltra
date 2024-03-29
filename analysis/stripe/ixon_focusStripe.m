function focus= ixon_focusStripe(x,y,z,stripe)
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

[xx,yy]= meshgrid(x,y);

% Fourier Filtering
% Nfft = 2^10;
% sPSF = 1.3161;% Gaussian radius of PSF in real space in camera pixels
% qPSF = sqrt(1/(4*pi*sPSF^2));% Gaussian radius of ideal PSF in momentum space
% dX = x(2)-x(1);
% f_max = 1/dX;
% f = 1/2*linspace(-f_max,f_max,Nfft);    
% [fxx,fyy]=meshgrid(f,f);
% zf = fftshift(fft2(z,Nfft,Nfft));

% lpf = ones(length(fyy),length(fyy));
% hpf = ones(length(fyy),length(fyy));

% lpf = exp(-(fxx.^2+fyy.^2)/(2*qPSF.^2));
% hpf = 1-exp(-(fxx.^2+fyy.^2)/(2*.01));

% zf = zf.*hpf;
% zf = zf.*lpf;

% z = abs(ifft2(zf,length(y),length(x)));
    
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


x_coms = [];
y_coms = [];
z_stripes = zeros(length(y),length(x),length(nVec));
stripe_sum = [];

xBs=zeros(length(nVec),2);
yBs = zeros(length(nVec),2);

z_stripes={};
x_stripes ={};
y_stripes = {};
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
    zThis = edgetaper(zThis,PSF);    
    zThis = imgaussfilt(zThis,1);

    z_stripes{end+1} = zThis;
    x_stripes{end+1} =xThis;
    y_stripes{end+1} = yThis;    
    [xx2,yy2]=meshgrid(xThis,yThis);    
    stripe_sum(end+1) = sum(zThis,'all');
    x_coms(end+1) = sum(xx2.*zThis,'all')/sum(zThis,'all');
    y_coms(end+1) = sum(yy2.*zThis,'all')/sum(zThis,'all'); 
end

z0 = sum(stripe_sum,'all');

%% Compute Focusing Score for each stripe
scores = zeros(length(z_stripes),1);
% figure(2000);

for n=1:length(z_stripes)
    z_to_analyze = z_stripes{n};  
    x=x_stripes{n};
    y=y_stripes{n};     
    % [xx,yy]=meshgrid(x,y);
    % nme = 1025;
    % f = 1/2*linspace(-f_max,f_max,nme);    
    % [fa,fb]=meshgrid(f,f);
    % zf_stripe = fftshift(fft2(z_to_analyze,nme,nme));    
    % zf_psf = exp(-(fa.^2+fb.^2)/(2*qPSF.^2));    
    % bb(n)=abs(sum(zf_psf.*zf_stripe,'all'))/sum(zf_stripe,'all');

    mag = 80;
    pxsize = 16;
    latt_spacing = 0.532;
    site_per_px = pxsize/mag/latt_spacing;
    sc = 2;
    eff_site_per_px = site_per_px/sc;
    zeff = imresize(z_to_analyze,eff_site_per_px);
    zeff = zeff/sum(zeff,'all');


    [sharpnessScore, map]=MLVSharpnessMeasure(zeff);  

%     subplot(length(z_stripes),2,n*2-1)
%     imagesc(x,y,map); axis equal tight
% 
% 
%      subplot(length(z_stripes),2,n*2)
%     imagesc(x,y,z_to_analyze); axis equal tight
    % [sharpnessScore, map]=MLVSharpnessMeasure(z_to_analyze);  
    scores(n) = sharpnessScore;
    if stripe_sum(n)<(z0*0.08)
        scores(n)=NaN;
        continue;
    end    
end

%% Fit Focus to ycom
inds = ~isnan(scores);

if sum(inds)>3

    yyy = scores(inds);
    xxx = y_coms(inds)';
    nnn = stripe_sum(inds)';
    pp = polyfit(xxx,yyy,2);

    myfit=fittype('a*x.^2+b*x+c','independent','x',...
        'coefficients',{'a','b','c'});
    fitopt = fitoptions(myfit);
    fitopt.StartPoint = [pp(1) pp(2) pp(3)];
    fitopt.Weights = nnn;

    fout = fit(xxx,yyy,myfit,fitopt);

    pp = [fout.a fout.b fout.c];
    y0 = -pp(2)/(2*pp(1));

else
    yyy = scores(inds);

    y0 = mean(yyy);   
    pp = [0 0 y0];
end
% 
% y_coms
% scores

%%
focus = struct;
focus.scores = scores;
focus.x_coms = x_coms;
focus.y_coms = y_coms;
focus.sums = stripe_sum;
focus.y0 = y0;
focus.poly = pp;
end

