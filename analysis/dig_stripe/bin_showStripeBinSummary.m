function [hF] = bin_showStripeBinSummary(bindata,xVar,opts)

if nargin <3
    opts = struct;
end

if ~isfield(opts,'nCenter')
    opts.nCenter = 90;
end

if ~isfield(opts,'ControlVariable')
   opts.ControlVariable='qgm_plane_uwave_frequency_offset';
end


co =  [0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];


BS = [bindata(1).BinStripe(1)];
for kk=2:length(bindata)
    BS(kk)=bindata(kk).BinStripe(1);
end
P = [bindata.Params];
X = [P.(xVar)];

alpha = [BS.ModDepth];
N = [BS.Counts];

% goodInds = ones(1,length(N));
goodInds = logical([alpha>.75].*[N>0.5e6]);
goodInds = logical([N>0.5e6]);

goodInds = logical(ones(length(bindata),1));
% 
P = P(goodInds);
X = X(goodInds);
BS = BS(goodInds);

%% Extract Local and Global Phase
phin_n = mod(2*pi*(opts.nCenter./[BS.Lambda])-[BS.Phase],2*pi);
% phin_global = mod([BS.Phase]-pi/2,2*pi);

Navg=2;
phin_global = unwrapPhaseTime(X,phin_n,Navg);


%% Focusing Positions
S=[];
C=[];
XX=[];
II=[];
MYC=[];

max_score = zeros(length(bindata),3);
focus_center = zeros(length(bindata),3);
% goodInds = logical([N>0.5e6]);

 for kk=1:length(bindata)
     try
     Zb = bindata(kk).LatticeBin.Zbin;
     Zb(isnan(Zb))=0;
     Zb(Zb<500) = 0;
     if sum(Zb,'all')<5e5
         continue;
     end
     
    ss = [bindata(kk).BinStripe(1).Scores];
    ss(isnan(ss))=0;
    cs=[bindata(kk).BinStripe(1).Centers];
    xs = repmat(X(kk),[length(cs),1]);
    
    n = length(xs);
    n = min([6 n]);
    
    
    XX=[XX;xs(1:n)];    
    [~,iii]=sort(ss,'descend');    
    MYC = [MYC; co(1:n,:)];
    S = [S; ss(iii(1:n))];
    C=[C;cs(iii(1:n))];
    
    focus_center(kk,1) = cs(iii(1));
    max_score(kk,1) = ss(iii(1));
    
    
    focus_center(kk,2) = cs(iii(2));
    max_score(kk,2) = ss(iii(2)); 
    
        
    focus_center(kk,3) = cs(iii(3));
    max_score(kk,3) = ss(iii(3)); 
    
    
     catch ME
        keyboard 
     end
 end
 



%%
hF = figure(4000);
hF.Color='w';
hF.Position = [100 20 1500 900];
clf

t=uicontrol('style','text','string',opts.FigLabel,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left');
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];
    
% Local Phase
subplot(2,4,1);
plot(X,phin_n/(2*pi),'ko','markerfacecolor',co(2,:),'linewidth',1);
xlabel(xVar,'interpreter','none');
if isequal(xVar,'ExecutionDate')
    datetick x
end
ylabel(['phase at n=' num2str(opts.nCenter) ' (2\pi)']);
ylim([0 1]);
xlim([min(X) max(X)]);
set(gca,'YTick',[0 .25 .5 .75 1]);
grid on
set(gca,'box','on','linewidth',1,'fontsize',10);

% Unwrapped Phase
subplot(2,4,2);
plot(X,phin_global/(2*pi),'ko','markerfacecolor',co(1,:),'linewidth',1);
xlabel(xVar,'interpreter','none');
if isequal(xVar,'ExecutionDate')
    datetick x
end
ylabel(['unwrapped phase at n=' num2str(opts.nCenter) ' (2\pi)']);
hold on
grid on
set(gca,'box','on','linewidth',1,'fontsize',10);
yL=get(gca,'YLim');
if (yL(2)-yL(1))<1
   ylim([0 1]); 
end
xlim([min(X) max(X)]);


% Focus Position
subplot(2,4,3);
plot(X,focus_center(:,1),'o-','markerfacecolor',co(1,:),'linewidth',1,'markeredgecolor',co(1,:)*.5);
xlabel(xVar,'interpreter','none');
hold on
plot(X,focus_center(:,2),'o','markerfacecolor',co(2,:),'linewidth',1,'markeredgecolor',co(2,:)*.5);
plot(X,focus_center(:,3),'o','markerfacecolor',co(3,:),'linewidth',1,'markeredgecolor',co(3,:)*.5);

if isequal(xVar,'ExecutionDate')
    datetick x
end
ylabel('stripe position (sites)');
set(gca,'box','on','linewidth',1,'fontsize',10);
xlim([min(X) max(X)]);
grid on

p1=plot(X,[BS.FocusCenterFit],'k-','linewidth',1);

legend(p1,{'center fit'},'fontsize',6,'location','northwest');

% Focus Scores
subplot(2,4,4);
plot(X,max_score(:,1),'o','markerfacecolor',co(1,:),'linewidth',1,'markeredgecolor',co(1,:)*.5);
xlabel(xVar,'interpreter','none');
hold on
plot(X,max_score(:,2),'o','markerfacecolor',co(2,:),'linewidth',1,'markeredgecolor',co(2,:)*.5);
plot(X,max_score(:,3),'o','markerfacecolor',co(3,:),'linewidth',1,'markeredgecolor',co(3,:)*.5);

if isequal(xVar,'ExecutionDate')
    datetick x
end
ylabel('scores (arb)');
set(gca,'box','on','linewidth',1,'fontsize',10);
xlim([min(X) max(X)]);
grid on


% Control Variable
subplot(2,4,5);
plot(X,[P.(opts.ControlVariable)],'ko','markerfacecolor',co(1,:),'linewidth',1);
xlabel(xVar,'interpreter','none');
if isequal(xVar,'ExecutionDate')
    datetick x
end
ylabel(opts.ControlVariable,'interpreter','none');
set(gca,'box','on','linewidth',1,'fontsize',10);
xlim([min(X) max(X)]);
grid on

% Wavelength
subplot(2,4,6);
plot(X,[BS.Lambda],'ko','markerfacecolor',co(1,:),'linewidth',1);
xlabel(xVar,'interpreter','none');
if isequal(xVar,'ExecutionDate')
    datetick x
end
ylabel('wavelength \lambda (sites)');
set(gca,'box','on','linewidth',1,'fontsize',10);
xlim([min(X) max(X)]);
grid on

% Duty Cycle
subplot(2,4,7);
plot(X,[BS.Duty]*100,'ko','markerfacecolor',co(1,:),'linewidth',1);
xlabel(xVar,'interpreter','none');
if isequal(xVar,'ExecutionDate')
    datetick x
end
ylabel('duty cycle (%)');
set(gca,'box','on','linewidth',1,'fontsize',10);
ylim([0 100]);
set(gca,'YTick',100*[0 .25 .5 .75 1]);
grid on
xlim([min(X) max(X)]);

% Modulation Depth
subplot(2,4,8);
plot(X,100*[BS.ModDepth],'ko','markerfacecolor',co(1,:),'linewidth',1);
xlabel(xVar,'interpreter','none');
if isequal(xVar,'ExecutionDate')
    datetick x
end
ylabel('Mod Depth %');
ylim([0 100]);
xlim([min(X) max(X)]);
set(gca,'YTick',100*[0 .25 .5 .75 1]);
grid on
set(gca,'box','on','linewidth',1,'fontsize',10);


% % Rsquare
% subplot(2,4,8);
% plot(X,[BS.RSquareStripe],'ko','markerfacecolor',co(1,:),'linewidth',1);
% xlabel(xVar,'interpreter','none');
% if isequal(xVar,'ExecutionDate')
%     datetick x
% end
% ylabel('RSquareStripe','interpreter','none');
% set(gca,'box','on','linewidth',1,'fontsize',10);
% xlim([min(X) max(X)]);




end


