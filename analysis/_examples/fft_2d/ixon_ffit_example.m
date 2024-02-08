function ixon_ffit_example
% 2D FFTs can be a bit funny (at least to CF), also where do the 2*pi's go?

% This peice of code provides intuitiion about our 2D FFT code.
ixondata=struct;
x=1:512;
y=1:512;

ixondata(1).Z=zeros(length(y),length(x));
ixondata(1).X=x;
ixondata(1).Y=y;

[xx,yy]=meshgrid(x,y);


%% Plane wave
theta=30;    % Angle in degrees
Lmin=10;
L=15; % wavelength
A=10;  % Amplitude

plane_wave=@(x,y) A*sin(2*pi/L*(cosd(theta)*x+sind(theta)*y));
%%

ixondata.Z=plane_wave(xx,yy)+0*randn(length(y),length(x));
ixondata=ixon_computeFFT(ixondata);

hF=figure;
hF.Color='w';
hF.Position(3:4)=[1000 400];

subplot(121)
imagesc(x,y,ixondata.Z);
xlabel('position (px)');
ylabel('position (px)');

axis equal tight
set(gca,'ydir','normal');
colorbar
colormap(purplemap)
caxis([-6 6]);

subplot(122)
imagesc(ixondata.fft_F,ixondata.fft_F,abs(ixondata.fft_Z));
xlabel('frequency (1/px)');
ylabel('frequency (1/px)');

axis equal tight
set(gca,'ydir','normal');
colorbar
caxis([0 A/2]);
colormap(purplemap)


xlim([-1 1]/Lmin);
ylim([-1 1]/Lmin);


%%
end

