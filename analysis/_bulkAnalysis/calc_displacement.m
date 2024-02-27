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

%% Calculate the displacement
wXDT = 2*pi*33; %2*pi*Hz
T = 20e-9; %K
G = 2*pi*14; %2*pi*Hz

d = [-15:0.01:15]*1e-6;
w = 2*pi*[40,90,100,110,120,70];
B = zeros([length(w),length(d)]);
A = zeros([length(w),length(d)]);
amp = zeros([length(w),length(d)]);
for ww = 1:length(w)
    
    B(ww,:) = -aL^2*m*wXDT^2*qfit_real(T,G,w(ww)).*d./(hbar*w(ww));
    A(ww,:) = -aL^2*m*wXDT^2*qfit_imag(T,G,w(ww)).*d./(hbar*w(ww));
    amp(ww,:) = sqrt(B(ww,:).^2 + A(ww,:).^2);
end

%% Plot Results
clf(10);
f10 = figure(10);
for loop = 1:length(w)
    subplot(2,3,loop);
    plot(d/1e-6,B(loop,:)/1e-6, 'DisplayName', 'Asin(\phi)')
    hold on;
    plot(d/1e-6,A(loop,:)/1e-6, 'DisplayName', 'Acos(\phi)')
    plot(d/1e-6,amp(loop,:)/1e-6, 'DisplayName', 'A')
    xlabel('trap displacement ($\mu m$)','Interpreter','latex')
    ylabel('amplitude/phase response ($\mu m$)','Interpreter','latex')
    title(['$\omega_{trap} = 55 \mathrm{ Hz}, \omega_{drive} = \mathrm{ }$' num2str(w(loop)/(2*pi)) '$\mathrm{ Hz}$'],'Interpreter','latex')
    legend('Location','south');
    xlim([-15,15])
    ylim([-3 3])
    grid();
end







