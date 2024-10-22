%% Constants
h = 6.626e-34;
hbar = h/(2*pi);
kb = 1.380649e-23;
aL = 527e-9;
amu = 1.6e-27;
m = 40*amu;

%% Define the lookup tables for R and E
global Rvalues;
Rvalues_unscaled = table2array(readtable('Rvalues_unscaled_64_4Hz_200.csv'));
Rvalues = -1j*(aL/pi)*Rvalues_unscaled;

global energies;
energies_Hz = importdata('EnergyHz_64_4Hz_200.txt');
energies = h*(energies_Hz);

%% Calculate the expected Fsum

tHz = 563.4109123332288;
tK = tHz*h/kb;

m0 = hbar^2/(2*aL^2*h*tHz);

T = [1:0.01:5]*tK;

STB = besseli(1,2.*tK./T)./besseli(0,2.*tK./T);

%% Calculated the FSum from spectra

TT = [1:0.1:5]*tK;
GG = 2*pi*14; %2*pi*Hz
fsum = [];
for loop = 1:length(TT)
    fun = @(w) qfit_real(TT(loop),GG,w).*(aL^2)./hbar;
    
    fsum(loop) = m0*(2/pi)*integral(fun,0,inf);
end

%% Susceptibility at omega = 0 fsum

% TT = [1:0.1:5]*tK;
% GG = 2*pi*30; %2*pi*Hz
% fsum = [];
% fsum2 = [];
% for loop = 1:length(TT)
%     fun = @(w) qfit_real(TT(loop),GG,w).*(aL^2)./hbar;
%     fun2 = @(w) qfit_real(TT(loop),GG,w).*(aL^2)./hbar./w.^2;
%     fsum(loop) = m0*(2/pi)*integral(fun,0,inf);
%     fsum2(loop) = (2*pi*64.4)^2*m0*(2/pi)*integral(fun2,0,inf);
% end

TT = [3]*tK;
GG = 2*pi*[10:5:80]; %2*pi*Hz
fsum = [];
fsum2 = [];
for loop = 1:length(GG)
    fun = @(w) qfit_real(TT,GG(loop),w).*(aL^2)./hbar;
    fun2 = @(w) qfit_real(TT,GG(loop),w).*(aL^2)./hbar./w.^2;
    fsum(loop) = m*(2/pi)*integral(fun,0,inf);
    fsum2(loop) = (2*pi*64.4)^2*m*(2/pi)*integral(fun2,0,inf);
end

%% Plot

f1 = figure(111);
clf(f1);
% plot(T./tK,STB,'DisplayName', 'S_{TB}');
hold on;
plot(GG,fsum,'DisplayName','$m \frac{2}{\pi}\int_0^{\infty} d\omega Re[\sigma(\omega)]$')
plot(GG,fsum2,'DisplayName','$m \omega_{\mathrm{Trap}}^2 \frac{2}{\pi}\int_0^{\infty} d\omega \frac{Re[\sigma(\omega)]}{\omega^2}$')
% plot(TT./tK,fsum,'DisplayName','S_{xx}')
% plot(TT./tK,fsum2,'DisplayName','S2_{xx}')
% xlabel('T/t', FontSize=16);
xlabel('\Gamma (s^{-1})', FontSize=16);
ylabel('$\frac{m}{N}S_{XX}$','Interpreter','latex',FontSize=16);
title('T/t = 3', FontSize=16)
% title('\Gamma = 2\pi \times 30 Hz',FontSize=16)
box on;
legend(FontSize=12,Interpreter = 'Latex');
%%
ww = 2*pi*[0:0.5:80] ;
fun3 = @(w) qfit_real(3*tK,0,w).*(aL^2)./hbar./w.^2;
f2 = figure(222);
clf(f2);
hold on;
plot(ww,(2*pi*64.4)^3*m0*fun3(ww))
