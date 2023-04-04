function out = stripeFFT(z,opt)
% stripeFFT
%
% do a 2D FFT on some data that likely has a fluoresnce stripe pattern.
% This 2D analyis very fast and can be used as an initial guess for full
% fits, or to be used as a quick measurement.


%% 
if nargin==0
    opt=struct;
    opt.doDebug = 0;
end

%% OPtions

maskR = 2;

R=20;
%%

% Make the image odd number so that the FFT has a nice zero point
if ~mod(size(z,1),2)
    z= z(2:end,:);
end
if ~mod(size(z,2),2)
    z= z(:,2:end);
end

% Center point of image
yc0 = (size(z,1)+1)/2;
xc0 = (size(z,2)+1)/2;

%% Compute 2D FFT

% Computte the FFT
zf = fft2(z);
zf = fftshift(zf);

% Find the norm
zf_abs = abs(zf);

% Peak Distribituion
zf_0 = max(max(zf_abs));

% Get rid of data below a 1e-3 threshold
zf_abs(zf_abs<zf_0*1e-3)=0;

% sampling frequency of ff1 is 1 pixel
fs = 1; 

% y frequency vector
ly = size(zf_abs,1);
fy = (-ly/2:ly/2-1)/ly*fs;

% x frequency vector
lx = size(zf_abs,2);
fx = (-lx/2:lx/2-1)/lx*fs;

%% Restrict FFT to low momentum modes

% Reduce area of inspection
ix = xc0 +[-R:1:R];
iy = yc0 +[-R:1:R];

%  K vectors in restricted ROI
fy = fy(iy);
fx = fx(ix);

% Data in restricted ROI
zf_abs = zf_abs(iy,ix);
zf = zf(iy,ix);



%% Compute moments and calculates covariance matrix to find angle
% See the following wiki page on how to caluclate rotations
% https://en.wikipedia.org/wiki/Image_moment

% Center point of image
yc0 = (size(zf_abs,1)+1)/2;
xc0 = (size(zf_abs,2)+1)/2;

% Mask out the center so the gaussian mode doesn't affect data
ix=xc0+[-maskR:maskR];
iy=yc0+[-maskR:maskR];

zf_abs_masked = zf_abs;
zf_abs_masked(iy,ix) = 0;  

% Get rid of maximal point
zf_abs_masked(zf_abs_masked<max(max(zf_abs_masked))*.2)=0;

% Compute the moments
m00 = computeMoment(zf_abs_masked,0,0);
m10 = computeMoment(zf_abs_masked,1,0); % m10
m01 = computeMoment(zf_abs_masked,0,1);
m11 = computeMoment(zf_abs_masked,1,1);
m20 = computeMoment(zf_abs_masked,2,0);
m02 = computeMoment(zf_abs_masked,0,2);

% Find Centroid
x_bar = m10/m00;
y_bar = m01/m00;

% computer central moment
u20 = m20/m00 - x_bar.^2;
u02 = m02/m00 - y_bar.^2;
u11 = m11/m00 - x_bar*y_bar;

theta = 0.5*atan(2*u11/(u20-u02));
%% Analyze Rotated Data

% Rotate the image, fft, and abs of fft
zr = imrotate(z,theta*180/pi,'bilinear','crop');
zfr_abs = imrotate(zf_abs,theta*180/pi,'bilinear','crop');
zfr = imrotate(zf,theta*180/pi,'bilinear','crop');

% Find the 3 biggest peaks
yp = sum(zfr_abs,2);

% [vals,p]= findpeaks(yp,'Npeaks',3);
[~,p]=islocalmax(yp);
[p,ind]=sort(p,'descend');
p=ind(2:3);

% keyboard
% Momentum in pixel of the peak
k = range(p)/2;

% Wavelength of fringe pattern
lambda = 1/(range(fy(1:2))*k);

% Calculate phase
zfr_1d = sum(zfr,2);
phi = angle(zfr_1d(p(end)));

% Compute center points of rotated image
xc = computeMoment(zr,1,0)/computeMoment(zr,0,0);
yc = computeMoment(zr,0,1)/computeMoment(zr,0,0);

% Calculate gaussian radius from the second moment in rotated basis
s1 = sqrt(computeMoment(zr,2,0)/computeMoment(zr,0,0)-xc^2);
s2 = sqrt(computeMoment(zr,0,2)/computeMoment(zr,0,0)-yc^2);


% Compute center points of original image
xc = computeMoment(z,1,0)/computeMoment(z,0,0);
yc = computeMoment(z,0,1)/computeMoment(z,0,0);

% Calculate the gaussian amplitude which is assumed to be the peak mode
n = max(max(zf_abs));
a=1/(2*s1*s2);
A = n*a/pi;

%% Output

out = struct;
out.A = A;
out.xc = xc;
out.yc = yc;
out.s1 = s1;
out.s2 = s2;
out.theta = theta;
out.lambda = lambda;
out.phi = phi;

%% Debug
if opt.doDebug
   
    f=figure(913);
    clf
    
    subplot(221)
    imagesc(z)
    axis equal tight
    caxis([0 A]);
    xlabel('x')
    ylabel('y')
    
    hold on    
    tt=linspace(0,2*pi,100);

    xa = sqrt(s1*s2)*cos(tt);
    ya = sqrt(s1*s2)*sin(tt);


    
    l = sqrt(s1*s2)*[-1.5 1.5];
    plot(xc+xa,yc+ya,'r-','linewidth',2)


    plot(xc+l*cos(theta+pi/2),yc+l*sin(theta+pi/2),'r-','linewidth',2)

    set(gca,'YDir','normal');

   
    subplot(222)
    imagesc(fx,fy,zf_abs)
    axis equal tight
    caxis([0 n*.25]);
    xlabel('kx')
    ylabel('ky')
    hold on
    
    r=linspace(-max(fx),max(fx),2);
    plot(r*cos(theta+pi/2),r*sin(theta+pi/2),'r-','linewidth',1)
    
       set(gca,'YDir','normal');
       
           subplot(223)
    imagesc(fx,fy,zf_abs_masked)
    axis equal tight
    caxis([0 n*.25]);
    xlabel('kx')
    ylabel('ky')
    hold on
    
    r=linspace(-max(fx),max(fx),2);
    plot(r*cos(theta+pi/2),r*sin(theta+pi/2),'r-','linewidth',1)
    
       set(gca,'YDir','normal');
       

end

end

function M = computeMoment(Z,ix,iy)

x=1:size(Z,2);
y=1:size(Z,1);

[xx,yy]=meshgrid(x,y);

M = sum(sum(xx.^(ix).*yy.^(iy).*Z));

end

