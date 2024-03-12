%% Bulk Conductivity Analysis and Spectrum Fitting
%Author: Frank Corapi
%Date: February 22nd, 2024

%% Select runs
% runs = [
%     2024 02 21 08;
%     2024 02 21 09;
%     2024 02 21 10;
%     2024 02 21 11;
%     2024 02 21 12;
%     2024 02 21 13;
%     2024 02 21 14;
%     2024 02 21 15;
%     2024 02 21 17;
%     2024 02 21 18;
%     2024 02 21 19;
%     2024 02 21 20;
%     2024 02 21 21;
%     2024 02 21 22;
%     2024 02 21 23;
%     2024 02 21 24;
%     2024 02 21 25;
%     ];

% runs = [
%     2024 02 22 04;
%     2024 02 22 05;
%     2024 02 22 06;
%     2024 02 22 07;
% 
%     ];

% % 196.5 G
% runs = [
%     2024 02 23 06;
%     2024 02 23 07;
%     2024 02 23 08;
%     2024 02 23 09;
%     2024 02 23 10;
%     2024 02 23 11;
%     2024 02 23 12;
%     2024 02 23 13;
%     2024 02 23 14;
%     2024 02 23 16;
%     2024 02 23 17;
%     2024 02 23 18;
%     2024 02 23 19;
%     2024 02 23 20;
%     2024 02 23 21;
%     2024 02 23 22;
%     2024 02 23 23;
%     2024 02 23 24;
%     2024 02 23 25;
%     ];

% % 190 G
% runs = [
%     2024 02 26 04;
%     2024 02 26 05;
%     2024 02 26 06;
%     2024 02 26 07;
%     2024 02 26 08;
%     2024 02 26 09;
%     2024 02 26 10;
%     2024 02 26 11;
%     2024 02 26 12;
%     2024 02 26 14;
%     2024 02 26 15;
%     2024 02 26 16;
%     2024 02 26 17;
%     2024 02 26 18;
%     2024 02 26 19;
%     2024 02 26 20;
%     2024 02 26 21;
%     2024 02 26 22;
%     2024 02 26 23;
%     2024 02 26 24;
%     ];

% % 201 G
% runs = [
%     2024 02 27 05;
%     2024 02 27 06;
%     2024 02 27 07;
%     2024 02 27 08;
%     2024 02 27 09;
%     2024 02 27 10;
%     2024 02 27 11;
%     2024 02 27 13;
%     2024 02 27 14;
%     2024 02 27 15;
%     2024 02 27 16;
%     2024 02 27 17;
%     2024 02 27 18;
%     2024 02 27 19;
%     2024 02 27 20;
%     2024 02 27 21;
% 
%     ];

% % 200.4 G
% runs = [
% 
%     2024 03 01 13;
%     2024 03 01 14;
%     2024 03 01 15;
%     2024 03 01 16;
%     2024 03 01 17;
%     2024 03 01 18;
%     2024 03 01 19;
%     2024 03 01 21;
%     2024 03 01 22;
%     2024 03 01 23;
%     2024 03 01 24;
%     2024 03 01 25;
%     2024 03 01 26;
%     2024 03 01 27;
%     2024 03 01 28;
%     2024 03 01 29;
%     ];

% 199.4 G
% runs = [
%     2024 03 04 05;
%     2024 03 04 06;
%     2024 03 04 07;
%     2024 03 04 08;
%     2024 03 04 09;
%     2024 03 04 10;
%     2024 03 04 11;
%     2024 03 04 13;
%     2024 03 04 14;
%     2024 03 04 15;
%     2024 03 04 16;
%     2024 03 04 17;
%     2024 03 04 18;
%     2024 03 04 19;
%     2024 03 04 20;
%     2024 03 04 21;
%     ];

% runs = [
%     2024 03 05 06;
%     2024 03 05 07;
%     2024 03 05 08;
%     2024 03 05 09;
%     2024 03 05 10;
%     2024 03 05 11;
%     2024 03 05 12;
%     2024 03 05 14;
%     2024 03 05 15;
%     2024 03 05 16;
%     2024 03 05 17;
%     2024 03 05 18;
%     2024 03 05 19;
%     2024 03 05 20;
%     2024 03 05 21;
%     2024 03 05 22;
%     ];

% runs = [
%     2024 03 06 07;
%     2024 03 06 08;
%     2024 03 06 09;
%     2024 03 06 10;
%     2024 03 06 11;
%     2024 03 06 12;
%     2024 03 06 13;
%     2024 03 06 15;
%     2024 03 06 16;
%     2024 03 06 17;
%     2024 03 06 18;
%     2024 03 06 19;
%     2024 03 06 20;
%     2024 03 06 21;
%     2024 03 06 22;
%     2024 03 06 23;
%     ];

% 200.8 G
% runs = [
% %     2024 03 07 06;
%     2024 03 07 07;
%     2024 03 07 08;
%     2024 03 07 09;
%     2024 03 07 10;
%     2024 03 07 11;
%     2024 03 07 12;
%     2024 03 07 14;
%     2024 03 07 15;
%     2024 03 07 16;
%     2024 03 07 17;
%     2024 03 07 18;
%     2024 03 07 19;
%     2024 03 07 20;
%     2024 03 07 21;
%     2024 03 07 22;
%     
%     ];

%200.9 G
runs = [
    2024 03 08 05;
    2024 03 08 06;
    2024 03 08 07;
    2024 03 08 08;
    2024 03 08 09;
    2024 03 08 10;
    2024 03 08 11;
    2024 03 08 13;
    2024 03 08 14;
    2024 03 08 15;
    2024 03 08 16;
    2024 03 08 17;
    2024 03 08 18;
    2024 03 08 19;
    2024 03 08 20;
    2024 03 08 21;
    
    ];

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
omega_ord = [];
field = [];
cond_real = [];
cond_real_err = [];
cond_imag = [];
cond_imag_err = [];
cond_real_ord = [];
cond_real_err_ord = [];
cond_imag_ord = [];
cond_imag_err_ord = [];
XsAvg = [];
YsAvg = [];
XsAvg_err = [];
YsAvg_err = [];
XsAvgsq_err = [];
YsAvgsq_err = [];

%Define constants
h = 6.626e-34; %Js
hbar = h/(2*pi); %Js
um = 1e-6; %m/um
kb = 1.380649e-23; %J/K
aL = 527e-9; %m
w_XDT = 2*pi*42.5; %2*pi*Hz
amu = 1.660538921e-27; %kg
m = 39.964008*amu; %kg
a_0 = 5.29177210903e-11; %m
w_T_xdir = 2*pi*66.8; % s^-1
w_T_ydir = 2*pi*59.7; % s^-1

%Experimental parameters
T = 80e-9; %K
G = 2*pi*18; %2*pi*Hz


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
    field(nn) = files{nn}.Params(1).conductivity_FB_field_maybe_calibrated;
    XsAvg(nn) = mean(files{nn}.Xs);
    YsAvg(nn) = mean(files{nn}.Ys);
    
    %Calculate standard deviation of cloud widths
    XsAvg_err(nn) = std(files{nn}.Xs);
    YsAvg_err(nn) = std(files{nn}.Ys);
    XsAvgsq_err(nn) = 2*std(files{nn}.Xs)*XsAvg(nn);
    YsAvgsq_err(nn) = 2*std(files{nn}.Ys)*YsAvg(nn);
    
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

%Sort the expected responses and conductivities
[freq_ord,freq_I] = sort(freq);
Cexp_ord = Cexp(freq_I);
Sexp_ord = Sexp(freq_I);
Aexp_ord = Aexp(freq_I);

omega_ord = omega(freq_I);
cond_real_ord = cond_real(freq_I);
cond_real_err_ord = cond_real_err(freq_I);
cond_imag_ord = cond_imag(freq_I);
cond_imag_err_ord = cond_imag_err(freq_I);

%% Experimental Conditions and Interaction Energies

% Feshbach Resonance
B_0 = 202.15; %G %free space feshbach resonance
deltaFB = 6.910; %G %width of FR
a_bg = 166.978*a_0; %m %Background scattering length
% Bfield = 200.91+0.1238; %G %Feshbach field
Bfield = field(1); %G %Feshbach field

as = a_bg*(1-deltaFB/(Bfield-B_0)); %Scattering length

% Lattice Condtions
W4 = 1.3798e19; %This is the integral of the ground band Wannier function to the 4th power, used to calculate U
g = 4*pi*as*hbar^2/m;

Us = (g)*W4/h; %Hz
tunneling = 563.4109123332288; %Hz
filling = 0.07; %spin-up atoms per site

nU2t = 2*pi*filling*Us^2/tunneling;

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
xlim([min(freq) max(freq)])
ylabel('amplitude/phase response (\mum)')
legend('S = A cos(\phi)','C = A sin(\phi)','A_{SC} = (S^2+C^2)^{1/2}',...
 'S expected','C expected','A expected',...
    'Location','northeast','FontSize',6)
grid on

subplot(122)
pA = errorbar(freq,A,Aerr,'o','Color',[.1 .6 .2],'MarkerFaceColor',[.1 .6 .2]);
% pA = errorbar(freq,abs(A),Aerr,'o','Color',[.1 .6 .2],'MarkerFaceColor',[.1 .6 .2]); % for negative
hold on
pA_SC2 = errorbar(freq,A_SC,A_SCerr,'ko','MarkerFaceColor','k');
pAexp2 = plot(freq_ord,Aexp_ord,'k--');
xlabel('frequency (Hz)')
xlim([min(freq) max(freq)])
ylabel('amplitude response (\mum)')
legend('A','A_{SC} = (S^2+C^2)^{1/2}','A expected')
grid on

% % Plot phase
% f6 = figure(666); % \m/
% clf;
% f6.WindowStyle = 'docked';
% f6.Color = 'w';
% 
% pphi = errorbar(freq,phi/pi,phierr/pi,'o','Color',[.1 .6 .2],'MarkerFaceColor',[.1 .6 .2]);
% % pphi = errorbar(freq,phi,phierr,'go','MarkerFaceColor','g'); % radians
% hold on
% pphiSC = errorbar(freq,phiSC/pi,phiSCerr/pi,'ko','MarkerFaceColor','k');
% % pphiSC = errorbar(freq,phiSC,phiSCerr,'ko','MarkerFaceColor','k'); % radians
% xlabel('frequency (Hz)')
% xlim([min(freq) max(freq)])
% ylabel('phase / \pi')
% legend('\phi','\phi_{SC} = arctan(C/S)')

%% Plot cloud width and Boltzmann temperature estimate

% Convert width to temperature
temp_const = m*w_T_xdir*w_T_ydir/kb/1e-9*1e-12; %nK/um^2

% Create blank docked figure with white background
f7 = figure(777);
clf;
f7.WindowStyle = 'docked';
f7.Color = 'w';

subplot(221)
pXsAvg = errorbar(freq,XsAvg,XsAvg_err,'o','Color',[.1 .2 .6],'MarkerFaceColor',[.1 .2 .6]);
xlabel('frequency (Hz)')
xlim([min(freq) max(freq)])
ylabel('x sigma (\mum)')

subplot(222)
pYsAvg = errorbar(freq,YsAvg,YsAvg_err,'o','Color',[.6 .2 .1],'MarkerFaceColor',[.6 .2 .1]);
xlabel('frequency (Hz)')
xlim([min(freq) max(freq)])
ylabel('y sigma (\mum)')

subplot(223)
pXsAvgsq = errorbar(freq,temp_const*(XsAvg.^2),temp_const*XsAvgsq_err,'o','Color',[.1 .2 .6],'MarkerFaceColor',[.1 .2 .6]);
xlabel('frequency (Hz)')
xlim([min(freq) max(freq)])
ylabel('temperature (nK)')
title('Boltzmann temperature x direction')

subplot(224)
pYsAvgsq = errorbar(freq,temp_const*(YsAvg.^2),temp_const*YsAvgsq_err,'o','Color',[.6 .2 .1],'MarkerFaceColor',[.6 .2 .1]);
xlabel('frequency (Hz)')
xlim([min(freq) max(freq)])
ylabel('temperature (nK)')
title('Boltzmann temperature y direction')

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
fout_real = fit(omega_ord(1:(end-5))',cond_real_ord(1:(end-5))',myfit_real,opt_real);
fout_real_c = confint(fout_real,lvl);
A_real_unc = (fout_real_c(2,1)-fout_real_c(1,1))/2;
B_real_unc = (fout_real_c(2,2)-fout_real_c(1,2))/2;
C_real_unc = (fout_real_c(2,3)-fout_real_c(1,3))/2;

opt_imag = fitoptions(myfit_imag);
opt_imag.StartPoint = [3600 2*pi*20 2*pi*50];
fout_imag = fit(omega_ord(1:(end-5))',cond_imag_ord(1:(end-5))',myfit_imag,opt_imag);
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
energies_Hz = importdata('EnergyHz_64_4Hz_200.txt');

global energies;
energies = h*(energies_Hz);

omega_pk = (energies(2)-energies(1))/hbar;

m_eff = m*(2*pi*64.4/omega_pk)^2;

global Rvalues;
Rvalues_unscaled = table2array(readtable('Rvalues_unscaled_64_4Hz_200.csv'));
Rvalues = -1j*(aL/pi)*Rvalues_unscaled;

myqfit_real = fittype(@(TT,GG,ww) qfit_real(TT,GG,ww), 'independent',{'ww'},...
        'coefficients',{'TT','GG'});

lvl = 0.667;
qopt_real = fitoptions(myqfit_real);
qopt_real.Display = 'iter';
qopt_real.StartPoint = [90e-9 2*pi*25];

qfout_real = fit(omega_ord(1:(end-5))',cond_real_ord(1:(end-5))',myqfit_real,qopt_real);
qfout_real_c=confint(qfout_real,lvl);
qT_real_unc = (qfout_real_c(2,1)-qfout_real_c(1,1))/2;
qG_real_unc = (qfout_real_c(2,2)-qfout_real_c(1,2))/2;

myqfit_imag = fittype(@(TT,GG,ww) qfit_imag(TT,GG,ww), 'independent',{'ww'},...
        'coefficients',{'TT','GG'});

qopt_imag = fitoptions(myqfit_imag);
qopt_imag.Display = 'iter';
qopt_imag.StartPoint = [90e-9 2*pi*25];


qfout_imag = fit(omega_ord(1:(end-5))',cond_imag_ord(1:(end-5))',myqfit_imag,qopt_imag);
qfout_imag_c=confint(qfout_imag,lvl);
qT_imag_unc = (qfout_imag_c(2,1)-qfout_imag_c(1,1))/2;
qG_imag_unc = (qfout_imag_c(2,2)-qfout_imag_c(1,2))/2;

%% Plot conductivity fits
ff = 0:0.1:250;
f4 = figure(444);
clf
f4.WindowStyle ='docked';
f4.Color = 'w';

subplot(121)
errorbar(freq_ord(1:(end-5)),cond_real_ord(1:(end-5)),cond_real_err_ord(1:(end-5)),'ko','markerfacecolor','k');
hold on
errorbar(freq_ord((end-4):end),cond_real_ord((end-4):end),cond_real_err_ord((end-4):end),...
'o','Color',[.1 .6 .2],'markerfacecolor',[.1 .6 .2]);

plot(ff,myfunc_real(fout_real.A,fout_real.B,fout_real.C,2*pi*ff))
plot(ff,qfit_real(qfout_imag.TT,qfout_imag.GG,2*pi*ff),'b--')
plot(ff,qfit_real(qfout_real.TT,qfout_real.GG,2*pi*ff),'color','b')
hold on
text(1,14,'$y = A\frac{\omega^2B}{(\omega^2-C^2)^2+w^2B^2}$', 'Interpreter','latex','color','r')
text(1,13,['$\frac{\Gamma}{2\pi} = \frac{B}{2\pi} = $' num2str(round(Gamma_real/(2*pi),2)) '$\pm$' num2str(round(B_real_unc/(2*pi),2)) ' Hz'], 'Interpreter','latex','color','r')
text(1,12,['$m^* = \frac{\hbar}{a_L^2A} = $ ' num2str(round(m_real/amu,2)) '$\pm$' num2str(round(m_real_unc/amu,2)) ' amu'], 'Interpreter','latex','color','r')
text(1,11,['$C/2\pi = $' num2str(round(fout_real.C/(2*pi),2)) '$\pm$' num2str(round(C_real_unc/(2*pi),2)) ' Hz'], 'Interpreter','latex','color', 'r')

text(90,14,['$\omega_{\mathrm{pk}} = 2\pi\times$' num2str(round(omega_pk/(2*pi),2)) ' Hz'], 'Interpreter','latex','color', 'b')
text(90,13,['$m^* = $ ' num2str(round(m_eff/amu,2)) ' amu'], 'Interpreter','latex','color','b')
text(90,12,['$\frac{\Gamma}{2\pi} = $' num2str(round(qfout_real.GG/(2*pi),2)) '$\pm$' num2str(round(qG_real_unc/(2*pi),2)) ' Hz'], 'Interpreter','latex','color','b')
text(90,11,['$T = $' num2str(round(qfout_real.TT/(1e-9),2)) '$\pm$' num2str(round(qT_real_unc/(1e-9),2)) ' nK'], 'Interpreter','latex','color','b')

xlabel('frequency (Hz)')
ylabel('real conductivity (\sigma/\sigma_0)')
title(['U = ' num2str(round(Us,2)) ' Hz' ', t = ' num2str(round(tunneling,2)) ' Hz']);
xlim([0 150])
ylim([-2 15])
grid on;

subplot(122)
errorbar(freq_ord(1:(end-5)),cond_imag_ord(1:(end-5)),cond_imag_err_ord(1:(end-5)),'ko','markerfacecolor','k');
hold on
errorbar(freq_ord((end-4):end),cond_imag_ord((end-4):end),cond_imag_err_ord((end-4):end),...
'o','Color',[.1 .6 .2],'markerfacecolor',[.1 .6 .2]);

plot(ff,myfunc_imag(fout_imag.A,fout_imag.B,fout_imag.C,2*pi*ff))
plot(ff,qfit_imag(qfout_real.TT,qfout_real.GG,2*pi*ff),'b--')
plot(ff,qfit_imag(qfout_imag.TT,qfout_imag.GG,2*pi*ff),'color','b')
hold on;
text(1,9,'$y = A\frac{\omega(\omega^2-C^2)}{(\omega^2-C^2)^2+w^2B^2}$', 'Interpreter','latex','color','r')
text(1,8,['$\frac{\Gamma}{2\pi} = \frac{B}{2\pi} = $' num2str(round(Gamma_imag/(2*pi),2)) '$\pm$' num2str(round(B_imag_unc/(2*pi),2)) ' Hz'], 'Interpreter','latex','color','r')
text(1,7,['$m^* = \frac{\hbar}{a_L^2A} = $ ' num2str(round(m_imag/amu,2)) '$\pm$' num2str(round(m_imag_unc/amu,2)) ' amu'], 'Interpreter','latex','color', 'r')
text(1,6,['$\frac{C}{2\pi} = $ ' num2str(round(fout_imag.C/(2*pi),2)) '$\pm$' num2str(round(C_imag_unc/(2*pi),2)) ' Hz'], 'Interpreter','latex','color', 'r')

text(90,9,['$\omega_{\mathrm{pk}} = 2\pi\times$' num2str(round(omega_pk/(2*pi),2)) ' Hz'], 'Interpreter','latex','color', 'b')
text(90,8,['$m^* = $ ' num2str(round(m_eff/amu,2)) ' amu'], 'Interpreter','latex','color', 'b')
text(90,7,['$\frac{\Gamma}{2\pi} = $' num2str(round(qfout_imag.GG/(2*pi),2)) '$\pm$' num2str(round(qG_imag_unc/(2*pi),2)) ' Hz'], 'Interpreter','latex','color', 'b')
text(90,6,['$T = $' num2str(round(qfout_imag.TT/(1e-9),2)) '$\pm$' num2str(round(qT_imag_unc/(1e-9),2)) ' nK'], 'Interpreter','latex','color', 'b')

xlabel('frequency (Hz)')
ylabel('imag conductivity (\sigma/\sigma_0)')

title(['U = ' num2str(round(Us,2)) ' Hz' ', t = ' num2str(round(tunneling,2)) ' Hz']);
xlim([0 150])
ylim([-5 10])
grid on;


