%% Experimental Conditions and Interaction Energies

% Feshbach Resonance
B_0 = 202.15; %G %free space feshbach resonance
deltaFB = 6.910; %G %width of FR
a_bg = 166.978*a_0; %m %Background scattering length
Bfield_list = [190,196.5,198.5,199.4,200,200.4,200.8,200.9,201]+0.1238; %G %Feshbach field

% Lattice Condtions
W4 = 1.3798e19; %This is the integral of the ground band Wannier function to the 4th power, used to calculate U
tunneling = 563.4109123332288; %Hz
filling = 0.07; %spin-up atoms per site
nU2t_list = zeros(numel(Bfield_list),1);

for bb=1:numel(Bfield_list)
    
    as = a_bg*(1-deltaFB/(Bfield_list(bb)-B_0)); %Scattering length
    g = 4*pi*as*hbar^2/m;
    Us = (g)*W4/h; %Hz
    nU2t_list(bb) = 2*pi*filling*Us^2/tunneling;
    
end
%% Plot conductivity summary

% xxxx = [290.6,593.6,1026.7,1545.5,2279.9,3229.6,5160.5,7010.3]; 
yyyy = 4*pi.*[8,13.59,9.77 13.51,14.77,18.48,19.77,23.41,18.41];
yyyy_err = 4*pi.*[0.88,1.11,1.25,1.83,2.37,3.53,5.01,6.87,2.95];

% constant temperature fit values
% yyyy = 4*pi.*[8.49,13.38,9.77,13.77,15.59,20,19.77,18.41];
% yyyy_err = 4*pi.*[1.55,1.6,1.25,2.39,3.01,3.92,5.01,4.55];

f8 = figure(888);
clf;
f8.WindowStyle = 'Docked';
f8.Color = 'w';
% errorbar(xxxx,yyyy,yyyy_err,'ko','MarkerFaceColor','k');
errorbar(nU2t_list,yyyy,yyyy_err,'ko','MarkerFaceColor','k');
xlabel('$nU^2/t\: (\mathrm{s}^{-1})$', 'Interpreter','latex')
ylabel('FWHM Dissipation Rate (s^{-1})')
% title('Constant temperature conductivity fits')