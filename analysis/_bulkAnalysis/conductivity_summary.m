%% Experimental Conditions and Interaction Energies

% Feshbach Resonance
B_0 = 202.15; %G %free space feshbach resonance
deltaFB = 6.910; %G %width of FR
a_bg = 166.978*a_0; %m %Background scattering length
Bfield_list = [190,196.5,198.5,199.4,200,200.4,200.8,200.9,201,201.1]+0.1238; %G %Feshbach field

% Lattice Condtions
W4 = 1.3798e19; %This is the integral of the ground band Wannier function to the 4th power, used to calculate U
tunneling = 563.4109123332288; %Hz
filling = 0.07; %spin-up atoms per site
Us_list = zeros(numel(Bfield_list),1);
nU2t_list = zeros(numel(Bfield_list),1);

for bb=1:numel(Bfield_list)
    
    as = a_bg*(1-deltaFB/(Bfield_list(bb)-B_0)); %Scattering length
    g = 4*pi*as*hbar^2/m;
    
    Us_list(bb) = (g)*W4/h; %Hz
    nU2t_list(bb) = 2*pi*filling*Us_list(bb)^2/tunneling;
    
end
%% Plot conductivity summary

% xxxx = [290.6,593.6,1026.7,1545.5,2279.9,3229.6,5160.5,7010.3]; 
% yyyy = 4*pi.*[8,13.59,9.77 13.51,14.77,18.48,19.77,18.41];
% yyyy_err = 4*pi.*[0.88,1.11,1.25,1.83,2.37,3.53,5.01,2.95];

% constant temperature fit values
real_2gamma = 4*pi.*[8.49,13.38,9.77,13.77,15.59,20,19.77,23.41,18.41,16.41];
real_2gamma_err = 4*pi.*[1.55,1.6,1.25,2.39,3.01,3.92,5.01,6.87,4.55,2.83];

imag_2gamma = 4*pi.*[14.78,20.33,12.78,12.6,10.02,15.23,10.88,17.17,15.23,22.77];
imag_2gamma_err = 4*pi.*[3.05,3.77,1.76,3.54,2.65,7.18,4.93,4.62,6.99,3.77];

f8 = figure(888);
clf;
f8.WindowStyle = 'Docked';
f8.Color = 'w';
% errorbar(xxxx,yyyy,yyyy_err,'ko','MarkerFaceColor','k');
errorbar(nU2t_list,real_2gamma,real_2gamma_err,'ko','MarkerFaceColor','k');
hold on
% errorbar(nU2t_list,imag_2gamma,imag_2gamma_err,'ro','MarkerFaceColor','r');
tau_ramp = 50;
yline(2e3/tau_ramp,'k--')
xlabel('$nU^2/t\: (\mathrm{s}^{-1})$', 'Interpreter','latex')
ylabel('dissipation rate (s^{-1})');
ylim([0 400]);
legend('$2\Gamma$ of Re$[\sigma]$','$2/\tau_{\mathrm{ramp}}$','interpreter','latex','fontsize',14)

% text(8000,2*pi*1e3/tau_ramp,['2\pi/(' num2str(tau_ramp) ' ms)'],...
%     'verticalalignment','bottom','horizontalalignment','left',...
%     'fontsize',12);
% title('Constant temperature conductivity fits')

%% Plot interaction energies

f9 = figure(999);
clf;
f9.WindowStyle = 'Docked';
f9.Color = 'w';

pUs = plot(nU2t_list,Us_list/tunneling,'ko','MarkerFaceColor','k');
xlabel('$nU^2/t\: (\mathrm{s}^{-1})$', 'Interpreter','latex')
ylabel('U/t')
title('Interaction energies')

%% Plot temperatures