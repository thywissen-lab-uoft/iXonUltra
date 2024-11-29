function [hF,out] = dig_ac_conductivity(digdata,opts)




if nargin ==1
    opts = struct;
end

if ~isfield(opts,'Ndelta_bound')
    opts.Ndelta_bound =  1.75;
end

if ~isfield(opts,'RemoveBadData')
    opts.RemoveBadData =  true;
end

if ~isfield(opts,'NumSpins')
    opts.NumSpins = 2;
end

Natoms = [digdata.Natoms];
X = [digdata.Xc_um];
Xs = [digdata.Xs_um];
Y = [digdata.Yc_um];
Ys = [digdata.Ys_um];

P = [digdata.Params];
T = [P.conductivity_mod_time];
Tr = [P.conductivity_mod_ramp_time];


% Mark bad inds
Nmed=median(Natoms);
Nbar = mean(Natoms);
Nstd = std(Natoms);
delta = (Natoms-Nbar)/Nstd;

if opts.RemoveBadData
    bad_inds=[abs(delta)>opts.Ndelta_bound];
else
    bad_inds=[];
end


% [Natoms,bad_inds] = rmoutliers([digdata.Natoms]);
% T(bad_inds)=[];
% Tr(bad_inds)=[];
% Ys(bad_inds)=[];
% Y(bad_inds)=[];
% Xs(bad_inds)=[];
% X(bad_inds)=[];
% Natoms(bad_inds)=[];

Ttot = T + Tr;

%%

hF = figure;
hF.Color='w';
hF.Position = [100 100 800 400];

hF.Name = 'Digital AC Conductivity';

if isfield(opts,'FigLabel') && ~isempty(opts.FigLabel)
    tFig=uicontrol('style','text','string',opts.FigLabel,...
        'units','pixels','backgroundcolor',...
        'w','horizontalalignment','left');
    tFig.Position(4)=tFig.Extent(4);
    tFig.Position(3)=hF.Position(3);
    tFig.Position(1:2)=[5 hF.Position(4)-tFig.Position(4)];
end    


ax1=subplot(2,3,[1 2]);
co=get(gca,'colororder');
plot(Ttot(~bad_inds),X(~bad_inds),'o','markerfacecolor',co(1,:),...
    'linewidth',1,'markeredgecolor',co(1,:)*.3);
hold on
plot(Ttot(bad_inds),X(bad_inds),'o','markerfacecolor',co(2,:),...
    'linewidth',1,'markeredgecolor',co(2,:)*.3);


xlabel('total modulation time (ms)');
ylabel('x center (um)');

if opts.RemoveBadData
    str = ['ignoring dN>' num2str(opts.Ndelta_bound) '\sigma'];
    text(.99,.01,str,'units','normalized','verticalalignment','bottom',...
        'HorizontalAlignment','right','fontsize',8);

end

ax2=subplot(2,3,4);
plot(Ttot(~bad_inds),Xs(~bad_inds),'o','markerfacecolor',co(1,:),...
    'linewidth',1,'markeredgecolor',co(1,:)*.3);
hold on
plot(Ttot(bad_inds),Xs(bad_inds),'o','markerfacecolor',co(2,:),...
    'linewidth',1,'markeredgecolor',co(2,:)*.3);
xlabel('total modulation time (ms)');
ylabel('x sigma (um)');

ax3=subplot(2,3,5);
plot(Ttot(~bad_inds),Natoms(~bad_inds),'o','markerfacecolor',co(1,:),...
    'linewidth',1,'markeredgecolor','k');
hold on
plot(Ttot(bad_inds),Natoms(bad_inds),'o','markerfacecolor',co(2,:),...
    'linewidth',1,'markeredgecolor','k');
ylabel('atom number');
xlabel('total modulation time (ms)');

%% Initial Fit Guess

f = unique([digdata.Params.conductivity_mod_freq]);
if length(f)>1;warning('multiple oscillation frequencies detected');end
f = mean(f);
f = f*1e-3;

Ag = 0.5*(max(X)-min(X));
x0 = median(X);

phiVec = linspace(0,-2*pi,100);
phi_sse = zeros(length(phiVec),1);
for kk=1:length(phiVec)
    phi_sse(kk) = sum((-Ag*sin(2*pi*f*Ttot + phiVec(kk)) - X).^2);
end
[~,ii] = min(phi_sse);
phi = phiVec(ii);
%%

sinephasefit = fittype(@(A,phi,x0,t) -A*sin(2*pi*f*t+phi)+x0,...
    'independent','t','coefficients',{'A','phi','x0'});
sinephasefit_opt = fitoptions(sinephasefit);
sinephasefit_opt.StartPoint = [Ag phi x0];
strFit1 = '$-A\sin(2\pi f t+\phi) + x_0$';


sinephasefit = fittype(@(A,phi,x0,v0,a0,t) ...
    -A*sin(2*pi*f*t+phi) + x0 + v0*t + 0.5*a0*t.^2,...
    'independent','t','coefficients',{'A','phi','x0','v0','a0'});
sinephasefit_opt = fitoptions(sinephasefit);
sinephasefit_opt.StartPoint = [Ag phi x0 0 0];
strFit1 = '$-A\sin(2\pi f t+\phi) + x_0 + v_0t +0.5a_0t^2$';
fout_sine = fit(Ttot',X',sinephasefit,sinephasefit_opt);

cint = confint(fout_sine,0.667);
Aerr = (cint(2,1)-cint(1,1))*0.5;
phierr = (cint(2,2)-cint(1,2))*0.5;
x0err = (cint(2,3)-cint(1,3))*0.5;
v0err = (cint(2,4)-cint(1,4))*0.5;
a0err = (cint(2,5)-cint(1,5))*0.5;

tbl_f1={['A(' char(956) 'm)'], [num2str(fout_sine.A,'%.2f') char(177) num2str(Aerr,'%.2f')];
[char(966) '(rad)'], [num2str(fout_sine.phi,'%.2f') char(177) num2str(phierr,'%.2f')];
['x' char(8320) '(' char(956) 'm)'], [num2str(fout_sine.x0,'%.1f') char(177) num2str(x0err,'%.1f')];
['v' char(8320) '(' char(956) 'm/ms)'], [num2str(fout_sine.v0,'%.1e') char(177) num2str(v0err,'%.1e')];
['a' char(8320) '(' char(956) 'm/ms)' char(178)], [num2str(fout_sine.a0,'%.1e') char(177) num2str(a0err,'%.1e')];
};
% strFit2 = ['$-S\sin(2\pi f t)-C\cos(2\pi f t) +$' newline ...
    % '$x_0 + v_0t +0.5a_0t^2$'];

%sinecosinefit = fittype(@(S,C,x0,v0,a0,t) ...
%    -S*sin(2*pi*f*t)-C*cos(2*pi*f*t) + x0 + v0*t + 0.5*a0*t.^2,...
%    'independent','t','coefficients',{'S','C','x0','v0','a0'});
%sinecosinefit_opt = fitoptions(sinecosinefit);
%sinecosinefit_opt.StartPoint = [Ag*cos(fout_sine.phi) Ag*sin(fout_sine.phi) x0 0 0];
%strFit2 = ['$-S\sin(2\pi f t)-C\cos(2\pi f t) + x_0 + v_0t +0.5a_0t^2$'];


sinecosinefit = fittype(@(S,C,x0,v0,t) ...
    -S*sin(2*pi*f*t)-C*cos(2*pi*f*t) + x0 + v0*t,...
    'independent','t','coefficients',{'S','C','x0','v0'});
sinecosinefit_opt = fitoptions(sinecosinefit);
sinecosinefit_opt.StartPoint = [Ag*cos(fout_sine.phi) Ag*sin(fout_sine.phi) x0 0];
if opts.RemoveBadData
    sinecosinefit_opt.Weights = double(~bad_inds);
end
strFit2 = ['$-S\sin(2\pi f t)-C\cos(2\pi f t) + x_0 + v_0t$'];




fout_sinecosine = fit(Ttot',X',sinecosinefit,sinecosinefit_opt);
cint2 = confint(fout_sinecosine,0.667);
Serr2 = (cint2(2,1)-cint2(1,1))*0.5;
Cerr2 = (cint2(2,2)-cint2(1,2))*0.5;
x0err2 = (cint2(2,3)-cint2(1,3))*0.5;
v0err2 = (cint2(2,4)-cint2(1,4))*0.5;
%a0err2 = (cint2(2,5)-cint2(1,5))*0.5;

%% Make Output


out = struct;
out.SourceDirecotry = digdata.SourceDirectory;
out.T = Ttot;
out.X = X;
out.Xs = Xs;
out.Y = Y;
out.Ys = Ys;
out.Natoms = Natoms;
out.Params = [digdata.Params];
out.Flags = [digdata.Flags];
out.Units = [digdata.Units];
out.FitSinePhase = fout_sine;
out.FitSineCosine = fout_sinecosine;

out.A = fout_sine.A;
out.phi = fout_sine.phi;
out.S = fout_sinecosine.S;
out.C = fout_sinecosine.C;
out.x0 = fout_sinecosine.x0;
out.v0 = fout_sinecosine.v0;
%out.a0 = fout_sinecosine.a0;


out.Aerr = Aerr;
out.phierr = phierr;
out.Serr =Serr2;
out.Cerr = Cerr2;
out.x0err = x0err2;
out.v0err = v0err2;
%out.a0err = a0err2;

Ysreal = Ys(~bad_inds);
Xsreal = Xs(~bad_inds);
Natomsreal = Natoms(~bad_inds);

Abar = (2*pi*mean(Ysreal)*mean(Xsreal)); % average area in um^-2
n = 0.532^2*mean(Natomsreal)/Abar/opts.NumSpins;% factor of 0.5 is for two spins
nerr = 0.532^2*std(Natomsreal)/Abar/opts.NumSpins;



out.YsBar               = mean(Ysreal);
out.YsBarErr            = std(Ysreal);
out.XsBar               = mean(Xsreal);
out.XsBarErr            = std(Xsreal);
out.Natoms              = mean(Natomsreal);
out.NatomsErr              = std(Natomsreal);
out.SpatialMaxUpDensity = n;
out.SpatialMaxUpDensityErr = nerr;

%% Make Table

% tbl_f2={['S(' char(956) 'm)'], [num2str(fout_sinecosine.S,'%.2f') char(177) num2str(Serr2,'%.2f')];
% ['C(' char(956) 'm)'], [num2str(fout_sinecosine.C,'%.2f') char(177) num2str(Cerr2,'%.2f')];
% ['x' char(8320) '(' char(956) 'm)'], [num2str(fout_sinecosine.x0,'%.1f') char(177) num2str(x0err2,'%.1f')];
% ['v' char(8320) '(' char(956) 'm/ms)'], [num2str(fout_sinecosine.v0,'%.1e') char(177) num2str(v0err2,'%.1e')];
% ['a' char(8320) '(' char(956) 'm/ms)' char(178)], [num2str(fout_sinecosine.a0,'%.1e') char(177) num2str(a0err2,'%.1e')];
% };
% 
% 
% tbl_f3={
% ['A(' char(956) 'm)'], [num2str(fout_sine.A,'%.2f') ' ' char(177) ' ' num2str(Aerr,'%.2f')];
% [char(966) '(rad)'], [num2str(fout_sine.phi,'%.2f') ' ' char(177) ' ' num2str(phierr,'%.2f')];
% ['S(' char(956) 'm)'], [num2str(fout_sinecosine.S,'%.2f') ' ' char(177) ' ' num2str(Serr2,'%.2f')];
% ['C(' char(956) 'm)'], [num2str(fout_sinecosine.C,'%.2f') ' ' char(177) ' ' num2str(Cerr2,'%.2f')];
% ['x' char(8320) '(' char(956) 'm)'], [num2str(fout_sinecosine.x0,'%.1f') ' ' char(177) ' ' num2str(x0err2,'%.1f')];
% ['v' char(8320) '(' char(956) 'm/ms)'], [num2str(fout_sinecosine.v0,'%.1e') ' ' char(177) ' ' num2str(v0err2,'%.1e')];
% ['a' char(8320) '(' char(956) 'm/ms' char(178) ')'], [num2str(fout_sinecosine.a0,'%.1e') ' ' char(177) ' ' num2str(a0err2,'%.1e')];
% };


% tbl_f3={
% ['S(' char(956) 'm)'], [num2str(fout_sinecosine.S,'%.2f') ' ' char(177) ' ' num2str(Serr2,'%.2f')];
% ['C(' char(956) 'm)'], [num2str(fout_sinecosine.C,'%.2f') ' ' char(177) ' ' num2str(Cerr2,'%.2f')];
% ['x' char(8320) '(' char(956) 'm)'], [num2str(fout_sinecosine.x0,'%.1f') ' ' char(177) ' ' num2str(x0err2,'%.1f')];
% ['v' char(8320) '(' char(956) 'm/ms)'], [num2str(fout_sinecosine.v0,'%.1e') ' ' char(177) ' ' num2str(v0err2,'%.1e')];
% ['a' char(8320) '(' char(956) 'm/ms' char(178) ')'], [num2str(fout_sinecosine.a0,'%.1e') ' ' char(177) ' ' num2str(a0err2,'%.1e')];
% ['N'], [num2str(round(mean(Natomsreal))) ' ' char(177) ' ' num2str(round(std(Natomsreal))) ];
% ['n' char(8320) 'up'], [num2str(n,'%.2f')  ' ' char(177) ' ' num2str(nerr,'%.2f')  ];
% [char(963) 'x'  '(' char(956) 'm)'], [num2str(out.XsBar,'%.2f')  ' ' char(177) ' ' num2str(out.XsBarErr,'%.2f')  ];
% [char(963) 'y'  '(' char(956) 'm)'], [num2str(out.YsBar,'%.2f')  ' ' char(177) ' ' num2str(out.YsBarErr,'%.2f')  ];
% };

tbl_f3={
['S(' char(956) 'm)'], [num2str(fout_sinecosine.S,'%.2f') ' ' char(177) ' ' num2str(Serr2,'%.2f')];
['C(' char(956) 'm)'], [num2str(fout_sinecosine.C,'%.2f') ' ' char(177) ' ' num2str(Cerr2,'%.2f')];
['x' char(8320) '(' char(956) 'm)'], [num2str(fout_sinecosine.x0,'%.1f') ' ' char(177) ' ' num2str(x0err2,'%.1f')];
['v' char(8320) '(' char(956) 'm/ms)'], [num2str(fout_sinecosine.v0,'%.1e') ' ' char(177) ' ' num2str(v0err2,'%.1e')];
['N'], [num2str(round(mean(Natomsreal))) ' ' char(177) ' ' num2str(round(std(Natomsreal))) ];
['n' char(8320) 'up'], [num2str(n,'%.2f')  ' ' char(177) ' ' num2str(nerr,'%.2f')  ];
[char(963) 'x'  '(' char(956) 'm)'], [num2str(out.XsBar,'%.2f')  ' ' char(177) ' ' num2str(out.XsBarErr,'%.2f')  ];
[char(963) 'y'  '(' char(956) 'm)'], [num2str(out.YsBar,'%.2f')  ' ' char(177) ' ' num2str(out.YsBarErr,'%.2f')  ];
};
tt=linspace(min(Ttot),max(Ttot),1e3);
%% Plot Fits

axes(ax1)
hold on
% pF1 = plot(tt,feval(fout_sine,tt),'r-','linewidth',2);
% pF2 = plot(tt,feval(fout_sinecosine,tt),'g--','linewidth',2);
% legend([pF1 pF2],{strFit1,strFit2},'interpreter','latex','fontsize',6,'location','best');

% pF1 = plot(tt,feval(fout_sine,tt),'r-','linewidth',2);
pF2 = plot(tt,feval(fout_sinecosine,tt),'g-','linewidth',2);
legend([pF2],{strFit2},'interpreter','latex','fontsize',7,'location','best');


ax4=subplot(1,3,3);
p = ax4.Position;
delete(ax4);
myTable = uitable('parent',hF,'units','normalized','pOsition',p);
myTable.Data = tbl_f3;
set(myTable,'rowname',{},'columnname',{},'fontsize',8,'Columnwidth',{70 100});
myTable.Position(3:4) = myTable.Extent(3:4);

%% Output


end

