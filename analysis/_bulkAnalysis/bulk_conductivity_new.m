%% Bulk Conductivity Analysis and Spectrum Fitting
%Author: Frank Corapi
%Date: February 22nd, 2024

%% Select runs
runs = [
    2024 02 21 08;
    2024 02 21 09;
    2024 02 21 10;
    2024 02 21 11;
    2024 02 21 12;
    2024 02 21 13;
    2024 02 21 14;
    2024 02 21 15;
    2024 02 21 17;
    2024 02 21 18;
    2024 02 21 19;
    2024 02 21 20;
    2024 02 21 21;
    2024 02 21 22;
    2024 02 21 23;
    2024 02 21 24;
    2024 02 21 25;
    ];

% runs = [
%     2024 02 22 04;
%     2024 02 22 05;
%     2024 02 22 06;
%     2024 02 22 07;
% 
%     ];


%% Find files and parameters and calculate the conductivity
dir_list = ixon_findRunDirectory(runs);

%Create empty lists
files = {};
A = [];
Aerr = [];
phi = [];
phierr = [];
S = [];
Serr = [];
C = [];
Cerr = [];
A_SC = [];
A_SCerr = [];
phiSC = [];
phiSCerr = [];
Aexp = [];
Sexp = [];
Cexp = [];
amp = [];
F = [];
freq = [];
omega = [];
cond_real = [];
cond_real_err = [];
cond_imag = [];
cond_imag_err = [];

%Define constants
h = 6.626e-34; %Js
hbar = h/(2*pi); %Js
um = 1e-6; %m/um
kb = 1.380649e-23; %J/K
aL = 527e-9; %m
w_XDT = 2*pi*33; %2*pi*Hz
amu = 1.660538921e-27; %kg
m = 39.964008*amu; %kg

%Experimental parameters
T = 19.37e-9; %K
G = 2*pi*26.14; %2*pi*Hz


for nn=1:length(dir_list)
    ixon_auto_dir = 0;
    imgdir = dir_list{nn};
    
    %Load Files
    files{nn}=load([imgdir filesep 'Figures' filesep 'conductivity_data.mat']);

    %Extract Parameters
    A(nn) = files{nn}.A; %um
    Aerr(nn) = files{nn}.Aerr; %um
    phi(nn) = files{nn}.phi; % rad
    phierr(nn) = files{nn}.phierr; % rad
    S(nn) = files{nn}.S; %um
    Serr(nn) = files{nn}.Serr; %um
    C(nn) = files{nn}.C; %um
    Cerr(nn) = files{nn}.Cerr; %um
    amp(nn) = 3.55*files{nn}.Params(1).conductivity_ODT2_mod_amp; %in um; Calibration to um is 3.55um/V
    F(nn) = m*(w_XDT^2)*amp(nn)*um; %N
    freq(nn) = files{nn}.Params(1).conductivity_mod_freq; %Hz
    omega(nn) = 2*pi*freq(nn); %2*pi*Hz
    
    %Calculate amplitude and phase from sine cosine fit
    A_SC(nn) = sqrt(S(nn)^2+C(nn)^2); %um
    A_SCerr(nn) = sqrt(S(nn)^2*Serr(nn)^2+C(nn)^2*Cerr(nn)^2)/A_SC(nn); %um
    
    phiSC(nn) = atan2(C(nn),S(nn)); %um
    phiSCerr(nn) = sqrt(C(nn)^2*Serr(nn)^2+S(nn)^2*Cerr(nn)^2)/A_SC(nn); %um
    
    %Calculate the expected displacement
    Cexp(nn) = (-aL^2*m*w_XDT^2*qfit_real(T,G,omega(nn))*amp(nn)*um/(hbar*omega(nn)))/um; %um
    Sexp(nn) = (-aL^2*m*w_XDT^2*qfit_imag(T,G,omega(nn))*amp(nn)*um/(hbar*omega(nn)))/um; %um
    Aexp(nn) = sqrt(Cexp(nn)^2 + Sexp(nn)^2); %um

    %Calculate the conductivity
    cond_real(nn) = -hbar*omega(nn)*(C(nn)*um)/(F(nn)*aL^2);
    cond_real_err(nn) = -hbar*omega(nn)*(Cerr(nn)*um)/(F(nn)*aL^2);

    cond_imag(nn) = -hbar*omega(nn)*(S(nn)*um)/(F(nn)*aL^2);
    cond_imag_err(nn) = -hbar*omega(nn)*(Serr(nn)*um)/(F(nn)*aL^2);
    
end

%Sort the expected responses
[freq_ord,freq_I] = sort(freq);
Cexp_ord = Cexp(freq_I);
Sexp_ord = Sexp(freq_I);
Aexp_ord = Aexp(freq_I);

%% Plot dig conductivity analysis fit parameters

% Plot amplitude fit parameters
f5 = figure(555);
clf;
f5.WindowStyle = 'docked';
f5.Color = 'w';

subplot(121)
pS = errorbar(freq,S,Serr,'ro','MarkerFaceColor','r');
hold on
pC = errorbar(freq,C,Cerr,'bo','MarkerFaceColor','b');
pA_SC = errorbar(freq,A_SC,A_SCerr,'ko','MarkerFaceColor','k');
pSexp = plot(freq_ord,Sexp_ord,'r--');
pCexp = plot(freq_ord,Cexp_ord,'b--');
pAexp = plot(freq_ord,Aexp_ord,'k--');
xlabel('frequency (Hz)')
ylabel('amplitude/phase response (\mum)')
legend('S = A cos(\phi)','C = A sin(\phi)','A_{SC} = (S^2+C^2)^{1/2}',...
 'S expected','C expected','A expected',...
    'Location','northeast','FontSize',6)

subplot(122)
pA = errorbar(freq,A,Aerr,'o','Color',[.1 .6 .2],'MarkerFaceColor',[.1 .6 .2]);
% pA = errorbar(freq,abs(A),Aerr,'o','Color',[.1 .6 .2],'MarkerFaceColor',[.1 .6 .2]); % for negative
hold on
pA_SC2 = errorbar(freq,A_SC,A_SCerr,'ko','MarkerFaceColor','k');
pAexp2 = plot(freq_ord,Aexp_ord,'k--');
xlabel('frequency (Hz)')
ylabel('amplitude response (\mum)')
legend('A','A_{SC} = (S^2+C^2)^{1/2}','A expected')

f6 = figure(666); % \m/
clf;
f6.WindowStyle = 'docked';
f6.Color = 'w';

pphi = errorbar(freq,phi/pi,phierr/pi,'o','Color',[.1 .6 .2],'MarkerFaceColor',[.1 .6 .2]);
% pphi = errorbar(freq,phi,phierr,'go','MarkerFaceColor','g'); % radians
hold on
pphiSC = errorbar(freq,phiSC/pi,phiSCerr/pi,'ko','MarkerFaceColor','k');
% pphiSC = errorbar(freq,phiSC,phiSCerr,'ko','MarkerFaceColor','k'); % radians
xlabel('frequency (Hz)')
ylabel('phase / \pi')
legend('\phi','\phi_{SC} = arctan(C/S)')


%% LRC Model Fit
myfunc_real = @(A,B,C,w) A*((B.*w.^2)./((w.^2-C.^2).^2+(w.*B).^2));
myfunc_imag = @(A,B,C,w) A*((w.^2-C.^2).*w./((w.^2-C.^2).^2+(w*B).^2));

myfit_real = fittype(@(A,B,C,w) myfunc_real(A,B,C,w),'independent',{'w'},...
        'coefficients',{'A','B','C'});
myfit_imag = fittype(@(A,B,C,w) myfunc_imag(A,B,C,w),'independent',{'w'},...
        'coefficients',{'A','B','C'});

lvl = 0.667;

opt_real = fitoptions(myfit_real);
opt_real.StartPoint = [3600 2*pi*20 2*pi*50];
fout_real = fit(omega',cond_real',myfit_real,opt_real);
fout_real_c = confint(fout_real,lvl);
A_real_unc = (fout_real_c(2,1)-fout_real_c(1,1))/2;
B_real_unc = (fout_real_c(2,2)-fout_real_c(1,2))/2;
C_real_unc = (fout_real_c(2,3)-fout_real_c(1,3))/2;

opt_imag = fitoptions(myfit_imag);
opt_imag.StartPoint = [3600 2*pi*20 2*pi*50];
fout_imag = fit(omega',cond_imag',myfit_imag,opt_imag);
fout_imag_c = confint(fout_real,lvl);
A_imag_unc = (fout_imag_c(2,1)-fout_imag_c(1,1))/2;
B_imag_unc = (fout_imag_c(2,2)-fout_imag_c(1,2))/2;
C_imag_unc = (fout_imag_c(2,3)-fout_imag_c(1,3))/2;

m_real = hbar/(fout_real.A*aL^2);
m_real_unc = m_real*A_real_unc/fout_real.A;
Gamma_real = fout_real.B;

m_imag = hbar/(fout_imag.A*aL^2);
m_imag_unc = m_imag*A_imag_unc/fout_imag.A;
Gamma_imag = fout_imag.B;

%% TDPT Fit (Full Quantum Model)
%See fit function definition at the end of the script

%Import energies and Rmn values from lookup tables
energies_Hz = importdata('EnergyHz_65.txt');

global energies;
energies = h*(energies_Hz);

omega_pk = (energies(2)-energies(1))/hbar;

m_eff = m*(2*pi*64.4/omega_pk)^2;

global Rvalues;
Rvalues_unscaled = table2array(readtable('Rvalues_unscaled_64_4Hz.csv'));
Rvalues = -1j*(aL/pi)*Rvalues_unscaled;


myqfit_real = fittype(@(TT,GG,ww) qfit_real(TT,GG,ww), 'independent',{'ww'},...
        'coefficients',{'TT','GG'});

lvl = 0.667;
qopt_real = fitoptions(myqfit_real);
qopt_real.Display = 'iter';
qopt_real.StartPoint = [30e-9 2*pi*20];

qfout_real = fit(omega',cond_real',myqfit_real,qopt_real);
qfout_real_c=confint(qfout_real,lvl);
qT_real_unc = (qfout_real_c(2,1)-qfout_real_c(1,1))/2;
qG_real_unc = (qfout_real_c(2,2)-qfout_real_c(1,2))/2;

myqfit_imag = fittype(@(TT,GG,ww) qfit_imag(TT,GG,ww), 'independent',{'ww'},...
        'coefficients',{'TT','GG'});

qopt_imag = fitoptions(myqfit_imag);
qopt_imag.Display = 'iter';
qopt_imag.StartPoint = [30e-9 2*pi*20];


qfout_imag = fit(omega',cond_imag',myqfit_imag,qopt_imag);
qfout_imag_c=confint(qfout_imag,lvl);
qT_imag_unc = (qfout_imag_c(2,1)-qfout_imag_c(1,1))/2;
qG_imag_unc = (qfout_imag_c(2,2)-qfout_imag_c(1,2))/2;

%% Plot conductivity fits
ff = 0:0.001:250;
f4 = figure(444);
clf
f4.WindowStyle ='docked';

subplot(121)
errorbar(freq,cond_real,cond_real_err,'ko','markerfacecolor','k');
hold on
plot(ff,myfunc_real(fout_real.A,fout_real.B,fout_real.C,2*pi*ff))
plot(ff,qfit_real(qfout_imag.TT,qfout_imag.GG,2*pi*ff),'b--')
plot(ff,qfit_real(qfout_real.TT,qfout_real.GG,2*pi*ff),'color','b')
hold on
text(1,3.3*2*pi,'$y = A\frac{\omega^2B}{(\omega^2-C^2)^2+w^2B^2}$', 'Interpreter','latex','color','r')
text(1,3*2*pi,['$\frac{\Gamma}{2\pi} = \frac{B}{2\pi} = $' num2str(round(Gamma_real/(2*pi),2)) '$\pm$' num2str(round(B_real_unc/(2*pi),2)) ' Hz'], 'Interpreter','latex','color','r')
text(1,2.7*2*pi,['$m^* = \frac{\hbar}{a_L^2A} = $ ' num2str(round(m_real/amu,2)) '$\pm$' num2str(round(m_real_unc/amu,2)) ' amu'], 'Interpreter','latex','color','r')
text(1,2.4*2*pi,['$C/2\pi = $' num2str(round(fout_real.C/(2*pi),2)) '$\pm$' num2str(round(C_real_unc/(2*pi),2)) ' Hz'], 'Interpreter','latex','color', 'r')

text(90,3.3*2*pi,['$\omega_{\mathrm{pk}} = 2\pi\times$' num2str(round(omega_pk/(2*pi),2)) ' Hz'], 'Interpreter','latex','color', 'b')
text(90,3*2*pi,['$m^* = $ ' num2str(round(m_eff/amu,2)) ' amu'], 'Interpreter','latex','color','b')
text(90,2.7*2*pi,['$\Gamma = 2\pi\times$' num2str(round(qfout_real.GG/(2*pi),2)) '$\pm$' num2str(round(qG_real_unc/(2*pi),2)) ' Hz'], 'Interpreter','latex','color','b')
text(90,2.4*2*pi,['$T = $' num2str(round(qfout_real.TT/(1e-9),2)) '$\pm$' num2str(round(qT_real_unc/(1e-9),2)) ' nK'], 'Interpreter','latex','color','b')

xlabel('frequency (Hz)')
ylabel('real conductivity (\sigma/\sigma_0)')
xlim([0 150])
ylim([-1.2 3.5]*2*pi)
grid on;

subplot(122)
errorbar(freq,cond_imag,cond_imag_err,'ko','markerfacecolor','k');
hold on;
plot(ff,myfunc_imag(fout_imag.A,fout_imag.B,fout_imag.C,2*pi*ff))
plot(ff,qfit_imag(qfout_real.TT,qfout_real.GG,2*pi*ff),'b--')
plot(ff,qfit_imag(qfout_imag.TT,qfout_imag.GG,2*pi*ff),'color','b')
hold on;
text(1,24,'$y = A\frac{\omega(\omega^2-C^2)}{(\omega^2-C^2)^2+w^2B^2}$', 'Interpreter','latex','color','r')
text(1,22,['$\frac{\Gamma}{2\pi} = \frac{B}{2\pi} = $' num2str(round(Gamma_imag/(2*pi),2)) '$\pm$' num2str(round(B_imag_unc/(2*pi),2)) ' Hz'], 'Interpreter','latex','color','r')
text(1,20,['$m^* = \frac{\hbar}{a_L^2A} = $ ' num2str(round(m_imag/amu,2)) '$\pm$' num2str(round(m_imag_unc/amu,2)) ' amu'], 'Interpreter','latex','color', 'r')
text(1,18,['$\frac{C}{2\pi} = $ ' num2str(round(fout_imag.C/(2*pi),2)) '$\pm$' num2str(round(C_imag_unc/(2*pi),2)) ' Hz'], 'Interpreter','latex','color', 'r')

text(90,24,['$\omega_{\mathrm{pk}} = 2\pi\times$' num2str(round(omega_pk/(2*pi),2)) ' Hz'], 'Interpreter','latex','color', 'b')
text(90,22,['$m^* = $ ' num2str(round(m_eff/amu,2)) ' amu'], 'Interpreter','latex','color', 'b')
text(90,20,['$\Gamma = 2\pi\times$' num2str(round(qfout_imag.GG/(2*pi),2)) '$\pm$' num2str(round(qG_imag_unc/(2*pi),2)) ' Hz'], 'Interpreter','latex','color', 'b')
text(90,18,['$T = $' num2str(round(qfout_imag.TT/(1e-9),2)) '$\pm$' num2str(round(qT_imag_unc/(1e-9),2)) ' nK'], 'Interpreter','latex','color', 'b')

xlabel('frequency (Hz)')
ylabel('imag conductivity (\sigma/\sigma_0)')
xlim([0 150])
ylim([-10 25])
grid on;


