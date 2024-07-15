
%% Constants
h=6.62607015e-34;
 hbar = 6.62607015e-34/(2*pi);

lam = (1054e-9*1054e-9*1064e-9)^(1/3);
k = 2*pi/lam;
m = 39.96399848*1.66053906660e-27;  

a_L = lam/2;
kL = k/2;

ERkHz = hbar*k^2/(2*m*2*pi)*1e-3;

% Harmonic trap parameters
w_XDT = 2*pi*42.5;
w_T_xdir = 2*pi*66.8; % s^-1
w_T_ydir = 2*pi*59.7; % s^-1

%% Load data

d = dig_radial_data;

N = d.AverageOccupation;
R = d.RadialVector;

% N = 2*N;

V = 1/2*m*w_XDT^2*a_L^2*R.^2;

% Parameters for 5ER lattice at 195G
U = h*1.563295893035468*1e3;
t = h*0.29369941733413896*1e3;

% Find max occupation

[Nmax,Imax] = max(N);
Rmax = R(Imax);
Vmax = V(Imax);

mu_0 = U/2 + Vmax;

% Calculate local chemical potential
mu = mu_0 - V;

% Fit to gaussian
gaussfit = fittype(@(N,sigma,x) ...
    N*exp(-(x-0.5).^2/(2*sigma^2)),...
    'independent','x','coefficients',{'N','sigma'});
gaussfit_opt = fitoptions(gaussfit);
gaussfit_opt.StartPoint = [Nmax 0.5];
fout_gauss = fit(mu/U,N,gaussfit,gaussfit_opt);


%% Plotting
clf

fmu = figure(100);
set(fmu,'color','w','Position',[50 50 1000 450]);

plot(mu/U,N,'o');

xlabel('$\mu/U$','Interpreter','latex');
    ylabel('Occupation ');
hold on;

x = -10:0.01:2;

plot(x,feval(fout_gauss,x));

hold off

%% Plotting


fmuV = figure(101);
set(fmuV,'color','w','Position',[50 50 1000 450]);

plot(mu/U,N.*(1-N),'o');

xlabel('$\mu/U$','Interpreter','latex');
ylabel('Variance');

figure(16)
plot(R,N,'o')







