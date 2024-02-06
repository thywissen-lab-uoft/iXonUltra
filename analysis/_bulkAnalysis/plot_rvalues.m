%% Plot the Rvalues in a 2D heatmap
Rvalues_unscaled = table2array(readtable('Rvalues_unscaled_64_4Hz.csv'));

start = 0;
final = 21;
n = [start:final];
m = [start:final];

f5 = figure(1001);
hm = heatmap(n,m,abs(Rvalues_unscaled((start+1):(final+1),(start+1):(final+1)))*(2/pi));
hm.YDisplayData=flip(hm.YDisplayData);
clim([0,23])
colormap(bone)
xlabel('m')
ylabel('n')
title('|R_{mn}| (a_L/2)')

%% Plot the real conductivity for various temperatures
freq = [0:200];

f10 = figure(1002);
plot(freq,qfit_real(300e-9,2*pi*16.6,2*pi*freq), 'DisplayName', 'T = 300 nK');
hold on;
plot(freq,qfit_real(70e-9,2*pi*16.6,2*pi*freq), 'DisplayName', 'T = 70 nK');
plot(freq,qfit_real(35e-9,2*pi*16.6,2*pi*freq), 'DisplayName', 'T = 35 nK');
plot(freq,qfit_real(20e-9,2*pi*16.6,2*pi*freq), 'DisplayName', 'T = 20 nK');
plot(freq,qfit_real(10e-9,2*pi*16.6,2*pi*freq), 'DisplayName', 'T = 10 nK');
xlabel('frequency (Hz)')
ylabel('real conductivity (\sigma/\sigma_0)')
title('$\Gamma/2\pi =$ 15 Hz, $\omega_0 = 55$ Hz, $V_L = 2.5 E_\mathrm{R}$','Interpreter','latex')

hold off;
legend;

%% Plot the imaginary cond, for various temperatures
freq = [0:0.1:200];

f1004 = figure(1004);
plot(freq,qfit_imag(300e-9,2*pi*16.6,2*pi*freq), 'DisplayName', 'T = 300 nK');
hold on;
plot(freq,qfit_imag(70e-9,2*pi*15,2*pi*freq), 'DisplayName', 'T = 70 nK');
plot(freq,qfit_imag(35e-9,2*pi*15,2*pi*freq), 'DisplayName', 'T = 35 nK');
plot(freq,qfit_imag(20e-9,2*pi*15,2*pi*freq), 'DisplayName', 'T = 20 nK');
plot(freq,qfit_imag(10e-9,2*pi*15,2*pi*freq), 'DisplayName', 'T = 10 nK');
xlabel('frequency (Hz)')
ylabel('imag conductivity (\sigma/\sigma_0)')
title('$\Gamma/2\pi =$ 15 Hz, $\omega_0 = 55$ Hz, $V_L = 2.5 E_\mathrm{R}$','Interpreter','latex')

hold off;
legend;

%% Find the peak frequency of the real conductivity
freq = [0:0.005:200];
T = [10:1:300]*1e-9;
fmax = zeros(1,length(T));
index = 1;
for loop = T
    q = qfit_real(loop,2*pi*15,2*pi*freq);
    [M,I] = max(q);
    fmax(index)=freq(I);
    index=index+1;
end
clf(f11);
f11 = figure(1003);
plot(T/1e-9,fmax)
xlabel('T (nK)')
ylabel('$\omega_{\mathrm{max}}/2\pi$ (Hz)','Interpreter','latex')
title('$\Gamma/2\pi =$ 15 Hz, $\omega_0 = 55$ Hz, $V_L = 2.5 E_\mathrm{R}$','Interpreter','latex')

%% Find the zero crossing of the real conductivity 

freq2 = [0:0.01:200];
T2 = [10:1:300]*1e-9;
fzero = zeros(1,length(T2));
index = 1;
for loop = T2
    q2 = qfit_imag(loop,2*pi*15,2*pi*freq2);
    zci = find(diff(sign(q2)));
    fzero(index)=freq2(zci(2));
    index=index+1;
end

%% Find the res frequency of the real conductivity
freqres = [0:0.01:200];
Tres = [10:1:300]*1e-9;
fres = zeros(1,length(Tres));
index = 1;
for loop = Tres
    q3 = qfit_imag(loop,2*pi*15,2*pi*freqres);
    [Mmax,Imax] = findpeaks(q3);
    [Mmin,Imin] = min(q3);
    I2 = round((Imax(1)+Imin)/2);
    fres(index)=freqres(I2);
    index=index+1;
end

%% Find the res frequency of the real conductivity
freqgrad = [0:0.01:200];
Tgrad = [10:1:300]*1e-9;
fgrad = zeros(1,length(Tgrad));
index = 1;
for loop = Tgrad
    q4 = qfit_imag(loop,2*pi*15,2*pi*freqgrad);
    [Mmax,Imax] = max(gradient(q4));
    fgrad(index)=freqgrad(Imax);
    index=index+1;
end

%% Plot
f1005 = figure(1005);
plot(T2/1e-9,fzero,'DisplayName', 'Imag Cond Zero')
hold on;
plot(T/1e-9,fmax,'DisplayName', 'Real Cond Peak')
plot(Tres/1e-9,fres,'DisplayName', 'Imag Cond Center')
plot(Tgrad/1e-9,fgrad,'DisplayName', 'Imag Cond Peak Slope')
hold off;
legend;
xlabel('T (nK)')
ylabel('$\omega_{\mathrm{res}}/2\pi$ (Hz)','Interpreter','latex')
title('$\Gamma/2\pi =$ 15 Hz, $\omega_0 = 55$ Hz, $V_L = 2.5 E_\mathrm{R}$','Interpreter','latex')
