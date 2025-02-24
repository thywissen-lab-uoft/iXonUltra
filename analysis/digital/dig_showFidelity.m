function [out, hF] = dig_showFidelity(digdata,opts)


if nargin ==1
    opts = struct;
end

if ~isfield(opts,'FigureNumber')
    opts.FigureNumber = 4002;
end

n1 = digdata.n1;
n2 = digdata.n2;

Z = zeros(numel(n2),numel(n1),2);

%% Get the fidelity data
for nn=1:size(digdata.Zdig,3)
    
    Z(:,:,1) = digdata.Zdig(:,:,nn);
    Z(:,:,2) = digdata.Zdig_img2(:,:,nn);
    
    out(nn) = dig_Fidelity(Z,n1,n2);

    
end

%% Initialize Graphics

hF = figure;
clf
hF.Color='w';
hF.Position= [100 100 1200 500];
hF.Name = 'Digital Fidelity';

if isfield(opts,'FigLabel') && ~isempty(opts.FigLabel)
    tFig=uicontrol('style','text','string',opts.FigLabel,...
        'units','pixels','backgroundcolor',...
        'w','horizontalalignment','left');
    tFig.Position(4)=tFig.Extent(4);
    tFig.Position(3)=hF.Position(3);
    tFig.Position(1:2)=[5 hF.Position(4)-tFig.Position(4)];
end    
co=get(gca,'colororder');

% keyboard

subplot(231)
plot(digdata.X,[out.Nlost],'o','markerfacecolor',co(1,:),...
'linewidth',1,'markeredgecolor',co(1,:)*.5);
ylabel('Number of lost atoms');
xlabel(digdata.xVar,'interpreter','none');
if isequal(digdata.xVar,'ExecutionDate')
    datetick x
end

subplot(232)
plot(digdata.X,[out.Nhop],'o','markerfacecolor',co(1,:),...
'linewidth',1,'markeredgecolor',co(1,:)*.5);
ylabel('Number of hopping events');
xlabel(digdata.xVar,'interpreter','none');
if isequal(digdata.xVar,'ExecutionDate')
    datetick x
end

subplot(234)
plot(digdata.X,[out.Nlost_percent].*100,'o','markerfacecolor',co(1,:),...
'linewidth',1,'markeredgecolor',co(1,:)*.5);
ylabel('Percent lost (%)');
xlabel(digdata.xVar,'interpreter','none');
if isequal(digdata.xVar,'ExecutionDate')
    datetick x
end

subplot(235)
plot(digdata.X,[out.Nhop_percent].*100,'o','markerfacecolor',co(1,:),...
'linewidth',1,'markeredgecolor',co(1,:)*.5);
ylabel('Percent hopped (%)');
xlabel(digdata.xVar,'interpreter','none');
if isequal(digdata.xVar,'ExecutionDate')
    datetick x
end

subplot(233)
plot(digdata.X,[out.N1],'o','markerfacecolor',co(1,:),...
'linewidth',1,'markeredgecolor',co(1,:)*.5);
hold on;
plot(digdata.X,[out.N2],'o','markerfacecolor',co(2,:),...
'linewidth',1,'markeredgecolor',co(2,:)*.5);
ylabel('Atom number');
xlabel(digdata.xVar,'interpreter','none');
if isequal(digdata.xVar,'ExecutionDate')
    datetick x
end



end