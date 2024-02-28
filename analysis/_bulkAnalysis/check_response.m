%% Constants
h = 6.626e-34;
hbar = h/(2*pi);
kb = 1.380649e-23;
aL = 527e-9;
amu = 1.660538921e-27;
m = 39.964008*amu;

%% Define the lookup tables for R and E
global Rvalues;
Rvalues_unscaled = table2array(readtable('Rvalues_unscaled_64_4Hz_200.csv'));
Rvalues = -1j*(aL/pi)*Rvalues_unscaled;

global energies;
energies_Hz = importdata('EnergyHz_64_4Hz_200.txt');
energies = h*(energies_Hz);

%% Trap parameters
wXDT = 2*pi*33; %2*pi*Hz
T = 20e-9; %K
G = 2*pi*80; %2*pi*Hz

%% Create data points

% x = [20,30,40,50,60,70,80,90,100,110,120];
% y = [0.45,0.38,0.30,0.29,0.40,0.61,0.89,1.22,1.60,2.02,2.51];

amp_desired = 0.65; %um
f = [20:5:150];
d = [];
for ff = 1:length(f)
    d(ff) = hbar*2*pi*f(ff)*amp_desired/(aL^2*m*wXDT^2*sqrt(qfit_real(T,G,2*pi*f(ff))^2+qfit_imag(T,G,2*pi*f(ff))^2))/3.55; %V
end

%% Fit the data
% myfunc = @(A,x0,w,x) A.*exp(-(x-x0).^2./w.^2) + 4;
myfunc = @(a,b,c,x) a.*(x-b).^2 + c;

myfit = fittype(@(a,b,c,x) myfunc(a,b,c,x),'independent',{'x'},...
        'coefficients',{'a','b','c'});

opt = fitoptions(myfit);
opt.StartPoint = [-1.5 20 40];
fout = fit(f',d',myfit,opt);

xx = 20:1:150;
yy = myfunc(fout.a,fout.b,fout.c,xx);


f33 = figure(33);
clf(f33);
plot(f,d,'ko')
hold on;
plot(xx,yy,'r')
% plot([20 40], [1.63,1.63],'r')
% plot([98 100], [4,4],'r')
xlabel('Drive Frequency (Hz)')
ylabel('Drive Amplitude (V)')
text(22,4.1,'y =a(x-b)^2 + c','color','r')
text(22,3.9,['a = ' num2str(fout.a)],'color','r')
text(22,3.7,['b = ' num2str(fout.b)],'color','r')
text(22,3.5,['c = ' num2str(fout.c)],'color','r')
