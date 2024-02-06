%% Define the lookup tables for R and E
global Rvalues;
Rvalues_unscaled = table2array(readtable('Rvalues_unscaled_55Hz.csv'));
Rvalues = -1j*(aL/pi)*Rvalues_unscaled;

global energies;
energies_Hz = importdata('EnergyHz_V55_2.txt');
energies = h*(energies_Hz);

%% Constants
h = 6.626e-34;
hbar = h/(2*pi);
kb = 1.380649e-23;
aL = 532e-9;

ff = 1:1:200;

creal = KKI2R(30e-9,2*pi*18,2*pi*ff);
cimag = KKR2I(30e-9,2*pi*18,2*pi*ff);


%% Plot
clf(f2);
f2 = figure(222);
plot(ff,creal,'DisplayName','Kramers-Kronig');
hold on;
plot(ff,qfit_real(30e-9,2*pi*18,2*pi*ff),'--','DisplayName','Calculated Real')
xlabel('frequency (Hz)')
ylabel('real conductivity (\sigma/\sigma_0)')
legend();

clf(f3);
f3 = figure(333);
plot(ff,cimag, 'DisplayName','Kramers-Kronig')
hold on;
plot(ff,qfit_imag(30e-9,2*pi*18,2*pi*ff),'--','DisplayName','Calculated Imag')
xlabel('frequency (Hz)')
ylabel('imag conductivity (\sigma/\sigma_0)')
legend();
%% Functions

function kk = KKI2R(T,G,w)
    
    eps = 1e-6;
    for loop = 1:length(w)
        fun = @(ww) qfit_imag(T,G,ww)./(ww-w(loop));
        
        kk1 = (1/pi)*integral(fun,-inf,w(loop)-eps);
        kk2 = (1/pi)*integral(fun,w(loop)+eps,inf);
    
        kk(loop) = kk1+kk2;
        disp(loop);
        
    end

end

function kk = KKR2I(T,G,w)
    
    eps = 1e-6;
    for loop = 1:length(w)
        fun2 = @(ww) qfit_real(T,G,ww)./(ww-w(loop));
        
        kk1 = (-1/pi)*integral(fun2,-inf,w(loop)-eps);
        kk2 = (-1/pi)*integral(fun2,w(loop)+eps,inf);
    
        kk(loop) = kk1+kk2;
        disp(loop)
    end

end