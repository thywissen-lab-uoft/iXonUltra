function [hF,out] = dig_bootstrap_ac_conductivity(digdata,opts)
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

%% Create figure

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

% if opts.RemoveBadData
%     str = ['ignoring dN>' num2str(opts.Ndelta_bound) '\sigma'];
%     text(.99,.01,str,'units','normalized','verticalalignment','bottom',...
%         'HorizontalAlignment','right','fontsize',8);
% end

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

keyboard

end

