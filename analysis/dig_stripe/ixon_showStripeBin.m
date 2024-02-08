function [] = ixon_showStripeBin(qgmdata,xVar,opts)

if nargin <3
    opts = struct;
end

if ~isfield(opts,'nCenter')
    opts.nCenter = 90;
end

if ~isfield(opts,'ControlVariable')
   opts.ControlVariable='qgm_plane_uwave_frequency_offset';
end

BS = [qgmdata.BinStripe];
P = [qgmdata.Params];
X = [P.(xVar)];

alpha = [BS.ModDepth];
N = [BS.Counts];

% goodInds = ones(1,length(N));
goodInds = logical([alpha>.75].*[N>0.5e6]);
goodInds = logical([N>0.5e6]);

goodInds = logical(ones(length(qgmdata),1));
% 
P = P(goodInds);
X = X(goodInds);
BS = BS(goodInds);

%% Extract Local and Global Phase
phin_n = mod(2*pi*(opts.nCenter./[BS.Lambda])-[BS.Phase],2*pi);
phin_global = mod([BS.Phase]-pi/2,2*pi);

%% Focusing Positions
S=[];
C=[];
XX=[];
II=[];
MYC=[];
co = get(gca,'colororder');

% goodInds = logical([N>0.5e6]);

 for kk=1:length(qgmdata)
    ss = [qgmdata(kk).BinStripe.Scores];
    ss(isnan(ss))=0;
    cs=[qgmdata(kk).BinStripe.Centers];
    xs = repmat(X(kk),[length(cs),1]);
    XX=[XX;xs(1:6)];    
    [~,iii]=sort(ss,'descend');    
    MYC = [MYC; co(1:6,:)];
    S = [S; ss(iii(1:6))];
    C=[C;cs(iii(1:6))];
 end


%%



hF = figure(4000);
hF.Color='w';
hF.Position = [100 100 900 900];
clf

subplot(4,2,1);
plot(X,phin_global/(2*pi),'ko','markerfacecolor',co(1,:),'linewidth',1);
hold on
plot(X,phin_n/(2*pi),'ko','markerfacecolor',co(2,:),'linewidth',1);
xlabel(xVar);
if isequal(xVar,'ExecutionDate')
    datetick x
end
ylabel('phase \phi (2\pi)');
ylim([0 1]);
xlim([min(X) max(X)]);
p=legend({'global',['\phi(site=' num2str(opts.nCenter) ')']},'orientation','horizontal',...
    'location','best','fontsize',8);
set(gca,'YTick',[0 .25 .5 .75 1]);
grid on
set(gca,'box','on','linewidth',1,'fontsize',10);

subplot(4,2,2);
plot(X,100*[BS.ModDepth],'ko','markerfacecolor',co(1,:),'linewidth',1);
xlabel(xVar);
if isequal(xVar,'ExecutionDate')
    datetick x
end
ylabel('Mod Depth %');
ylim([0 100]);
xlim([min(X) max(X)]);
set(gca,'YTick',100*[0 .25 .5 .75 1]);
grid on
set(gca,'box','on','linewidth',1,'fontsize',10);

subplot(4,2,3);
Sn = S/max(S);
Sn = Sn*50;
Sn(Sn==0)=1;
scatter(XX,C,Sn,MYC,'filled')
xlabel(xVar);
if isequal(xVar,'ExecutionDate')
    datetick x
end
ylabel('focus position');
set(gca,'box','on','linewidth',1,'fontsize',10);
xlim([min(X) max(X)]);

subplot(4,2,4);
plot(X,[BS.Counts],'ko','markerfacecolor',co(1,:),'linewidth',1);
xlabel(xVar);
if isequal(xVar,'ExecutionDate')
    datetick x
end
ylabel('counts');
yL=get(gca,'ylim');
set(gca,'Ylim',[0 yL(2)]);
set(gca,'box','on','linewidth',1,'fontsize',10);
xlim([min(X) max(X)]);

subplot(4,2,5);
plot(X,[BS.Lambda],'ko','markerfacecolor',co(1,:),'linewidth',1);
xlabel(xVar);
if isequal(xVar,'ExecutionDate')
    datetick x
end
ylabel('wavelength \lambda (sites)');
set(gca,'box','on','linewidth',1,'fontsize',10);
xlim([min(X) max(X)]);

subplot(4,2,6);
plot(X,[BS.Duty]*100,'ko','markerfacecolor',co(1,:),'linewidth',1);
xlabel(xVar);
if isequal(xVar,'ExecutionDate')
    datetick x
end
ylabel('duty cycle (%)');
set(gca,'box','on','linewidth',1,'fontsize',10);
ylim([0 100]);
set(gca,'YTick',100*[0 .25 .5 .75 1]);
grid on
xlim([min(X) max(X)]);


subplot(4,2,7);
plot(X,[P.(opts.ControlVariable)],'ko','markerfacecolor',co(1,:),'linewidth',1);
xlabel(xVar);
if isequal(xVar,'ExecutionDate')
    datetick x
end
ylabel(xVar);
set(gca,'box','on','linewidth',1,'fontsize',10);
xlim([min(X) max(X)]);



end


