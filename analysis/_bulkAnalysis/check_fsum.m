%% Constants
h = 6.626e-34;
hbar = h/(2*pi);
kb = 1.380649e-23;
aL = 527e-9;

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



%% Plot

f1 = figure(111);
clf(f1);
plot(T./tK,STB,'DisplayName', 'S_{TB}');
hold on;
plot(TT/tK,fsum,'DisplayName','S_{int}')
xlabel('T/t');
ylabel('$\frac{m_0^*}{N}S_{XX}^{TB}$','Interpreter','latex');
legend();
