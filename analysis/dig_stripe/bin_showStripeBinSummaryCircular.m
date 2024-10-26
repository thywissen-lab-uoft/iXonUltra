function hF = bin_showStripeBinSummaryCircular(bindata,xVar,opts)

if nargin <3
    opts = struct;
end

if ~isfield(opts,'nCenter')
    opts.nCenter = [110 110];
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

BS = [bindata(1).BinStripeCircular(1)];
for kk=2:length(bindata)
    BS(kk)=bindata(kk).BinStripeCircular(1);
    Phi(kk) = BS(kk).PhaseFunc(opts.nCenter(1),opts.nCenter(2));
end
P = [bindata.Params];
X = [P.(xVar)];

%% Get Values
Lambda=[BS.Lambda];
Theta = [BS.Theta];
Radius= [BS.Radius];

%% Fix Domain of phi to be be around 0 +-0.5

Phi=mod(Phi+pi,2*pi)-pi;

Navg=2;
PhiUnwrap = unwrapPhaseTime(X,Phi,Navg);

%% Smush Data


n1=bindata(1).LatticeBin(1).n1;

Zsum = zeros(length(n1),length(bindata));


for nn=1:length(BS)
    Zb=bindata(nn).LatticeBin(1).Zbin;
    Zb(isnan(Zb))=0;
    Zb(isinf(Zb))=0;

    Zbrot = imrotate(Zb,-Theta(nn),'crop');

    Zs = sum(Zbrot,1);
    Zs = Zs(:);
    Zs = Zs/sum(Zs);
    Zsum(:,nn)=Zs;
    

end


%%
hF = figure(4000);
hF.Color='w';
hF.Position = [100 20 1500 600];
clf

 ca = [1 1 1];
cb = [0.6 0 .5];
cc = [linspace(ca(1),cb(1),1000)' ...
    linspace(ca(2),cb(2),1000)' linspace(ca(3),cb(3),1000)'];
colormap(hF,cc);
                
                
t=uicontrol('style','text','string',opts.FigLabel,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left');
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];
    
% Local Phase
subplot(2,4,1);
plot(X,Phi/(2*pi),'ko','markerfacecolor',co(2,:),'linewidth',1);
xlabel(xVar,'interpreter','none');
if isequal(xVar,'ExecutionDate')
    datetick x
end
ylabel(['\phi(' num2str(opts.nCenter(1)) ',' num2str(opts.nCenter(2)) ') (2\pi)']);
ylim([-.5 .5]);
xlim([min(X) max(X)]);
set(gca,'YTick',[-.5 -.25 0 .25 .5]);
grid on
set(gca,'box','on','linewidth',1,'fontsize',10);

% Unwrapped Phase
subplot(2,4,2);
plot(X,PhiUnwrap/(2*pi),'ko','markerfacecolor',co(1,:),'linewidth',1);
xlabel(xVar,'interpreter','none');
if isequal(xVar,'ExecutionDate')
    datetick x
end
ylabel(['unwrap \phi(' num2str(opts.nCenter(1)) ',' num2str(opts.nCenter(2)) ') (2\pi)']);
hold on
grid on
set(gca,'box','on','linewidth',1,'fontsize',10);
yL=get(gca,'YLim');
% if (yL(2)-yL(1))<1
%    ylim([0 1]); 
% end
xlim([min(X) max(X)]);


% Image
subplot(2,4,3);
imagesc(X,n1,Zsum)
set(gca,'YDir','normal','box','on','linewidth',1);
hold on
xlabel(xVar,'interpreter','none');

if isequal(xVar,'ExecutionDate')
    datetick x
end
ylabel('position (site)');
set(gca,'box','on','linewidth',1,'fontsize',10);
xlim([min(X) max(X)]);
grid on


% Control Variable
subplot(2,4,4);
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
subplot(2,4,5);
plot(X,Lambda,'ko','markerfacecolor',co(1,:),'linewidth',1);
xlabel(xVar,'interpreter','none');
if isequal(xVar,'ExecutionDate')
    datetick x
end
ylabel('wavelength \lambda (sites)');
set(gca,'box','on','linewidth',1,'fontsize',10);
xlim([min(X) max(X)]);
grid on

% Rotation
subplot(2,4,6);
plot(X,Theta/pi*180,'ko','markerfacecolor',co(1,:),'linewidth',1);
xlabel(xVar,'interpreter','none');
if isequal(xVar,'ExecutionDate')
    datetick x
end
ylabel('rotation \theta (deg)');
set(gca,'box','on','linewidth',1,'fontsize',10);
xlim([min(X) max(X)]);
grid on


% Radius
subplot(2,4,7);
plot(X,Radius,'ko','markerfacecolor',co(1,:),'linewidth',1);
xlabel(xVar,'interpreter','none');
if isequal(xVar,'ExecutionDate')
    datetick x
end
ylabel('radius of curvature R (sites)');
set(gca,'box','on','linewidth',1,'fontsize',10);
xlim([min(X) max(X)]);
grid on


% Unwrapped Phase
subplot(2,4,8);
% plot(X,PhiUnwrap/(2*pi),'ko','markerfacecolor',co(1,:),'linewidth',1);
histogram(Phi/(2*pi),linspace(-1,1,25)*.5);
xlabel('\phi (2\pi)');
ylabel('occurences');
hold on
grid on
set(gca,'box','on','linewidth',1,'fontsize',10);
xlim([-.5 .5]);

end

