function [hF,hF2,out] = dig_bootstrap_ac_conductivity(digdata,opts)
%% Define missing arguments
if nargin ==1
    opts = struct;
end

if ~isfield(opts,'RemoveBadData')
    opts.RemoveBadData =  true;
end

if ~isfield(opts,'NumSpins')
    opts.NumSpins = 2;
end

if ~isfield(opts,'SaveFigs')
    opts.SaveFigs = 1;
end

%% Re-define variables

Natoms = [digdata.Natoms];
X = [digdata.Xc_um];
X = X';
Xs = [digdata.Xs_um];
Y = [digdata.Yc_um];
Ys = [digdata.Ys_um];

P = [digdata.Params];
T = [P.conductivity_mod_time];
Tr = [P.conductivity_mod_ramp_time];

Ttot = T + Tr;

%% Mark bad data

Nmed=median(Natoms);
Nbar = mean(Natoms);
Nstd = std(Natoms);

b1 = [Natoms>Nmed*1.5];
b2 = [Natoms<Nmed*0.5];
bad_inds = logical(b1)|logical(b2);

%% Read in the drive frequency

%Read the frequency of oscillation
f = unique([digdata.Params.conductivity_mod_freq]);
if length(f)>1;warning('multiple oscillation frequencies detected');end
f = mean(f);
f = f*1e-3;

%These guesses are not used
% Ag = 0.5*(max(X)-min(X));
% x0 = median(X);
% 
% phiVec = linspace(0,-2*pi,100);
% phi_sse = zeros(length(phiVec),1);
% for kk=1:length(phiVec)
% %     keyboard
%     phi_sse(kk) = sum((-Ag*sin(2*pi*f*Ttot + phiVec(kk)) - X).^2);
% %     phi_sse(kk) = sum((-Ag*sin(2*pi*f*Ttot + phiVec(kk)) - X.').^2);
% end
% 
% [~,ii] = min(phi_sse);
% phi_g = phiVec(ii); 

%% Post select bad data points based off atom number

data = [Ttot(~bad_inds)', X(~bad_inds)'];

%% Import normal least squares fit results
cond_data_temp = load([opts.saveDir filesep 'conductivity_data.mat']);

S = cond_data_temp.S;
C = cond_data_temp.C;
A = cond_data_temp.A;
phi = cond_data_temp.phi;
x0 = cond_data_temp.x0;
v0 = cond_data_temp.v0;

Serr = cond_data_temp.Serr;
Cerr = cond_data_temp.Cerr;
Aerr = cond_data_temp.Aerr;
phierr = cond_data_temp.phierr;
x0err = cond_data_temp.x0err;
v0err = cond_data_temp.v0err;

Natoms_mean = cond_data_temp.Natoms;
SpatialMaxUpDensity = cond_data_temp.SpatialMaxUpDensity;
XsBar = cond_data_temp.XsBar;
YsBar = cond_data_temp.YsBar;

Natoms_std = cond_data_temp.NatomsErr;
SpatialMaxUpDensityErr = cond_data_temp.SpatialMaxUpDensityErr;
XsBarErr = cond_data_temp.XsBarErr;
YsBarErr = cond_data_temp.YsBarErr;

%% Use fit phi as phase guess
phi_g = phi;

%% Bootstrap Fitting for SC model

qphi = opts.QPD_phi;
sinecosine = @(P,t) -P(1)*sin(2*pi*f*t+qphi)-P(2)*cos(2*pi*f*t+qphi) + P(3) + P(4)*t;

P_guess = [S C x0 v0];

nBootstraps = 1e4; %number of bootstrap iterations
[bootstat, bootsam] = bootstrp(nBootstraps, @fitModel_SC, data);

S_boot = bootstat(:,1);
C_boot = bootstat(:,2);
x0_boot = bootstat(:,3);
v0_boot = bootstat(:,4);

%% Fit the SC bootstrap data with a normal distribution
S_norm = fitdist(S_boot,'Normal');
C_norm = fitdist(C_boot,'Normal');
x0_norm = fitdist(x0_boot,'Normal');
v0_norm = fitdist(v0_boot,'Normal');

alpha = 1 - 0.95; %1 - confidence interval
S_ci = paramci(S_norm,'Alpha',alpha);
C_ci = paramci(C_norm,'Alpha',alpha);
x0_ci = paramci(x0_norm,'Alpha',alpha);
v0_ci = paramci(v0_norm,'Alpha',alpha);

S_mu = S_norm.mu;
C_mu = C_norm.mu;
x0_mu = x0_norm.mu;
v0_mu = v0_norm.mu;

S_sigma = S_norm.sigma;
C_sigma = C_norm.sigma;
x0_sigma = x0_norm.sigma;
v0_sigma = v0_norm.sigma;

%% Bootstrap Fitting for Aphi model

qphi = opts.QPD_phi;
sinephase = @(P,t) -P(1)*sin(2*pi*f*t+P(2)+qphi) + P(3) + P(4)*t;

P_Aphi_guess = [A phi_g x0 v0];

nBootstraps = 1e4; %number of bootstrap iterations
[bootstat2, bootsam2] = bootstrp(nBootstraps, @fitModel_Aphi, data);

A_boot = bootstat2(:,1);
phi_boot = bootstat2(:,2);
x02_boot = bootstat2(:,3);
v02_boot = bootstat2(:,4);

%% Fit the bootstrap data with a normal distribution
A_norm = fitdist(A_boot,'Normal');
phi_norm = fitdist(phi_boot,'Normal');
x02_norm = fitdist(x02_boot,'Normal');
v02_norm = fitdist(v02_boot,'Normal');

alpha = 1 - 0.95; %1 - confidence interval
A_ci = paramci(A_norm,'Alpha',alpha);
phi_ci = paramci(phi_norm,'Alpha',alpha);
x02_ci = paramci(x02_norm,'Alpha',alpha);
v02_ci = paramci(v02_norm,'Alpha',alpha);

A_mu = A_norm.mu;
phi_mu = phi_norm.mu;
x02_mu = x02_norm.mu;
v02_mu = v02_norm.mu;

A_sigma = A_norm.sigma;
phi_sigma = phi_norm.sigma;
x02_sigma = x02_norm.sigma;
v02_sigma = v02_norm.sigma;

%% Create table of data

tbl_f3={['Bootstrap SC: '],[''];
['S(' char(956) 'm)'], [num2str(S_mu,'%.2f') ' ' char(177) ' ' num2str(S_sigma,'%.2f')];
['C(' char(956) 'm)'], [num2str(C_mu,'%.2f') ' ' char(177) ' ' num2str(C_sigma,'%.2f')];
['x' char(8320) '(' char(956) 'm)'], [num2str(x0_mu,'%.1f') ' ' char(177) ' ' num2str(x0_sigma,'%.1f')];
['v' char(8320) '(' char(956) 'm/ms)'], [num2str(v0_mu,'%.1e') ' ' char(177) ' ' num2str(v0_sigma,'%.1e')];
['Bootstrap Aphi: '],[''];
['A(' char(956) 'm)'], [num2str(A_mu,'%.2f') ' ' char(177) ' ' num2str(A_sigma,'%.2f')];
[char(966) ' (rad.)'], [num2str(phi_mu,'%.2f') ' ' char(177) ' ' num2str(phi_sigma,'%.2f')];
['x' char(8320) '(' char(956) 'm)'], [num2str(x02_mu,'%.1f') ' ' char(177) ' ' num2str(x02_sigma,'%.1f')];
['v' char(8320) '(' char(956) 'm/ms)'], [num2str(v02_mu,'%.1e') ' ' char(177) ' ' num2str(v02_sigma,'%.1e')];
['Fit: '],[''];
['S(' char(956) 'm)'], [num2str(S,'%.2f') ' ' char(177) ' ' num2str(Serr,'%.2f')];
['C(' char(956) 'm)'], [num2str(C,'%.2f') ' ' char(177) ' ' num2str(Cerr,'%.2f')];
['x' char(8320) '(' char(956) 'm)'], [num2str(x0,'%.1f') ' ' char(177) ' ' num2str(x0err,'%.1f')];
['v' char(8320) '(' char(956) 'm/ms)'], [num2str(v0,'%.1e') ' ' char(177) ' ' num2str(v0err,'%.1e')];
['A(' char(956) 'm)'], [num2str(A,'%.2f') ' ' char(177) ' ' num2str(Aerr,'%.2f')];
[char(966) ' (rad.)'], [num2str(phi,'%.2f') ' ' char(177) ' ' num2str(phierr,'%.2f')];
['Density Info: '],[''];
['N'], [num2str(round(Natoms_mean)) ' ' char(177) ' ' num2str(round(Natoms_std)) ];
['n' char(8320) 'up'], [num2str(SpatialMaxUpDensity,'%.2f')  ' ' char(177) ' ' num2str(SpatialMaxUpDensityErr,'%.2f')  ];
[char(963) 'x'  '(' char(956) 'm)'], [num2str(XsBar,'%.2f')  ' ' char(177) ' ' num2str(XsBarErr,'%.2f')  ];
[char(963) 'y'  '(' char(956) 'm)'], [num2str(YsBar,'%.2f')  ' ' char(177) ' ' num2str(YsBarErr,'%.2f')  ];
};

%% Create figure showing bootstrap and fit results

hF = figure;
hF.Color='w';
hF.Position = [100 100 1400 600];

hF.Name = 'Digital AC Conductivity Bootstrap';

if isfield(opts,'FigLabel') && ~isempty(opts.FigLabel)
    tFig=uicontrol('style','text','string',opts.FigLabel,...
        'units','pixels','backgroundcolor',...
        'w','horizontalalignment','left');
    tFig.Position(4)=tFig.Extent(4);
    tFig.Position(3)=hF.Position(3);
    tFig.Position(1:2)=[5 hF.Position(4)-tFig.Position(4)];
end    

tt=linspace(min(Ttot),max(Ttot),1e3);
ax1=subplot(2,6,[1 2]);
co=get(gca,'colororder');
plot(Ttot(~bad_inds),X(~bad_inds),'o','markerfacecolor',co(1,:),...
    'linewidth',1,'markeredgecolor',co(1,:)*.3);
hold on
plot(Ttot(bad_inds),X(bad_inds),'o','markerfacecolor',co(2,:),...
    'linewidth',1,'markeredgecolor',co(2,:)*.3);
pF1 = plot(tt,sinecosine([S_mu,C_mu,x0_mu,v0_mu],tt), 'r','LineWidth',2);
pF2 = plot(tt,sinephase([A_mu,phi_mu,x02_mu,v02_mu],tt), 'k:','LineWidth',2);
pF3 = plot(tt,sinecosine([S,C,x0,v0],tt), 'g--','LineWidth',2);
strFit = {'Bootstrap SC', 'Bootstrap $A\phi$', '$-S\sin(2\pi f t)-C\cos(2\pi f t) + x_0 + v_0t$'};
legend([pF1 pF2 pF3],strFit,'interpreter','latex','fontsize',7,'location','best');

xlabel('total modulation time (ms)');
ylabel('x center (um)');

ax2=subplot(2,6,7);
plot(Ttot(~bad_inds),Xs(~bad_inds),'o','markerfacecolor',co(1,:),...
    'linewidth',1,'markeredgecolor',co(1,:)*.3);
hold on
plot(Ttot(bad_inds),Xs(bad_inds),'o','markerfacecolor',co(2,:),...
    'linewidth',1,'markeredgecolor',co(2,:)*.3);
xlabel('total modulation time (ms)');
ylabel('x sigma (um)');

ax3=subplot(2,6,8);
plot(Ttot(~bad_inds),Natoms(~bad_inds),'o','markerfacecolor',co(1,:),...
    'linewidth',1,'markeredgecolor','k');
hold on
plot(Ttot(bad_inds),Natoms(bad_inds),'o','markerfacecolor',co(2,:),...
    'linewidth',1,'markeredgecolor','k');
ylabel('atom number');
xlabel('total modulation time (ms)');

ax4 = subplot(2,6,3);
histfit(S_boot);
xlabel('S (\mum)')
ylabel('Occurences')
text(0.02,0.95,['\mu: ' num2str(round(S_mu,2)) '\pm' num2str(round(abs(S_ci(1,1) - S_ci(2,1))/2,2)) newline ...
    '\sigma: ' num2str(round(S_sigma,3)) '\pm' num2str(round(abs(S_ci(1,2) - S_ci(2,2))/2,4))],'Units','normalized', 'BackgroundColor','white');

ax5 = subplot(2,6,4);
histfit(C_boot);
xlabel('C (\mum)')
ylabel('Occurences')
text(0.02,0.95,['\mu: ' num2str(round(C_mu,2)) '\pm' num2str(round(abs(C_ci(1,1) - C_ci(2,1))/2,2)) newline ...
    '\sigma: ' num2str(round(C_sigma,3)) '\pm' num2str(round(abs(C_ci(1,2) - C_ci(2,2))/2,4))],'Units','normalized','BackgroundColor','white')

ax6 = subplot(2,6,9);
histfit(A_boot);
xlabel('A (\mum)')
ylabel('Occurences')
text(0.02,0.95,['\mu: ' num2str(round(A_mu,2)) '\pm' num2str(round(abs(A_ci(1,1) - A_ci(2,1))/2,2)) newline ...
    '\sigma: ' num2str(round(A_sigma,3)) '\pm' num2str(round(abs(A_ci(1,2) - A_ci(2,2))/2,4))],'Units','normalized', 'BackgroundColor','white');

ax7 = subplot(2,6,10);
histfit(phi_boot);
xlabel('\phi (rad.)')
ylabel('Occurences')
text(0.02,0.95,['\mu: ' num2str(round(phi_mu,2)) '\pm' num2str(round(abs(phi_ci(1,1) - phi_ci(2,1))/2,2)) newline ...
    '\sigma: ' num2str(round(phi_sigma,3)) '\pm' num2str(round(abs(phi_ci(1,2) - phi_ci(2,2))/2,4))],'Units','normalized', 'BackgroundColor','white');

ax8 = subplot(2,6,5);
histfit(x0_boot);
xlabel('x0\_SC (\mum)')
ylabel('Occurences')
text(0.02,0.95,['\mu: ' num2str(round(x0_mu,2)) '\pm' num2str(round(abs(x0_ci(1,1) - x0_ci(2,1))/2,2)) newline ...
    '\sigma: ' num2str(round(x0_sigma,3)) '\pm' num2str(round(abs(x0_ci(1,2) - x0_ci(2,2))/2,4))],'Units','normalized','BackgroundColor','white')

ax9 = subplot(2,6,11);
histfit(v0_boot);
xlabel('v0\_SC (\mum/ms)')
ylabel('Occurences')
text(0.02,0.95,['\mu: ' num2str(round(v0_mu,2)) '\pm' num2str(round(abs(v0_ci(1,1) - v0_ci(2,1))/2,2)) newline ...
    '\sigma: ' num2str(round(v0_sigma,3)) '\pm' num2str(round(abs(v0_ci(1,2) - v0_ci(2,2))/2,4))],'Units','normalized','BackgroundColor','white')

ax10=subplot(2,6,[6 12]);
p = ax10.Position;
delete(ax10);
myTable = uitable('parent',hF,'units','normalized','position',p);
myTable.Data = tbl_f3;
set(myTable,'rowname',{},'columnname',{},'fontsize',8,'Columnwidth',{90 100});
myTable.Position(3:4) = myTable.Extent(3:4);


%% Calculate correlations

SC_corr = mean((S_boot-S_mu).*(C_boot-C_mu))/(S_sigma*C_sigma);
Aphi_corr = mean((A_boot-A_mu).*(phi_boot-phi_mu))/(A_sigma*phi_sigma);


%% Create figure showing bootstrap correlations for SC and Aphi fits

hF2 = figure;
hF2.Color='w';
hF2.Position = [100 100 1400 600];

hF2.Name = 'Digital AC Conductivity Bootstrap Correlations';

if isfield(opts,'FigLabel') && ~isempty(opts.FigLabel)
    tFig=uicontrol('style','text','string',opts.FigLabel,...
        'units','pixels','backgroundcolor',...
        'w','horizontalalignment','left');
    tFig.Position(4)=tFig.Extent(4);
    tFig.Position(3)=hF2.Position(3);
    tFig.Position(1:2)=[5 hF2.Position(4)-tFig.Position(4)];
end 

ax10=subplot(1,2,1);
scatter(C_boot,S_boot,'k','MarkerFaceColor','k');
text(0.1,0.9,['\rho = ' num2str(round(SC_corr,4))],'Units','normalized','BackgroundColor','white','FontSize',12)
xlim([-3 3])
ylim([-3 3])
box on
grid on
xlabel('C (\mum)')
ylabel('S (\mum)')

ax20=subplot(1,2,2);
scatter(A_boot,phi_boot,'k','MarkerFaceColor','k');
text(0.1,0.9,['\rho = ' num2str(round(Aphi_corr,4))],'Units','normalized','BackgroundColor','white','FontSize',12)
% xlim([0 1.4])
ylim([-2*pi 2*pi])
box on
grid on
xlabel('A (\mum)')
ylabel('\phi (rad.)')



%% Make Output

out = struct;
out.SourceDirectory = digdata.SourceDirectory;
out.T = Ttot;
out.X = X;
out.Xs = Xs;
out.Y = Y;
out.Ys = Ys;
out.Natoms = Natoms_mean;
out.Params = [digdata.Params];
out.Flags = [digdata.Flags];
out.Units = [digdata.Units];

out.S = S_mu;
out.C = C_mu;
out.x0 = x0_mu;
out.v0 = v0_mu;
out.A = A_mu;
out.phi = phi_mu;
out.x02 = x02_mu;
out.v02 = v02_mu;

out.Serr = S_sigma;
out.Cerr = C_sigma;
out.x0err = x0_sigma;
out.v0err = v0_sigma;
out.A = A_sigma;
out.phi = phi_sigma;
out.x02 = x02_sigma;
out.v02 = v02_sigma;

out.SC_corr = SC_corr;
out.Aphi_corr = Aphi_corr;

out.YsBar               = YsBar;
out.YsBarErr            = YsBarErr ;
out.XsBar               = XsBar;
out.XsBarErr            = XsBarErr;
out.NatomsErr           = Natoms_std;
out.SpatialMaxUpDensity = SpatialMaxUpDensity;
out.SpatialMaxUpDensityErr = SpatialMaxUpDensityErr;


%% Define functions

function fitParams = fitModel_SC(data)
    options = optimset('Display','off');
    fitParams = lsqcurvefit(sinecosine, P_guess, data(:,1), data(:,2), [], [], options);
end

function fitParams = fitModel_Aphi(data)
    options = optimset('Display','off');
    fitParams = lsqcurvefit(sinephase, P_Aphi_guess, data(:,1), data(:,2), [0 (-2*pi+phi_g) -inf -inf], [20 (2*pi+phi_g) inf inf], options);
end


end

