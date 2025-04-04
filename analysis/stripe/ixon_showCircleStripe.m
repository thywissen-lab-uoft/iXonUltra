function hF = ixon_showCircleStripe(stripe,xVar,opts)

if nargin <3
    opts = struct;
end

if ~isfield(opts,'nCenter')
    opts.nCenter = [110 110];
end

if ~isfield(opts,'ControlVariable')
   opts.ControlVariable='qgm_plane_uwave_frequency_offset';
end

   opts.ControlVariable='f_offset';


co =  [0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

% stripe = [ixondata.StripeCircular];


for kk=1:length(stripe)
    Phi(kk) = stripe(kk).PhaseFunc(opts.nCenter(1),opts.nCenter(2));
end
P = [stripe.Params];
X = [P.(xVar)];

%% Get Values
Lambda=[stripe.Lambda];
Theta = [stripe.Theta];
Radius= [stripe.Radius];

%% Fix Domain of phi to be be around 0 +-0.5

Phi=mod(Phi+pi,2*pi)-pi;

Navg=3;
PhiUnwrap = unwrapPhaseTime(X,Phi,Navg);

%% Smush Data


% Xvec = ixondata(1).X;
% Yvec = ixondata(1).Y;

L = size(stripe(1).ZRotatedSum,1);


ZsAll = zeros(L,length(stripe));
for nn=1:length(stripe)
    % Z = sum(ixondata(nn).Z,3);
    % Z(isnan(Z))=0;
    % Z(isinf(Z))=0;
    % Zrot = imrotate(Z,-Theta(nn),'crop');
    % Zs = sum(Zrot,1);
    % Zs = Zs(:);
    % ZsAll(:,nn) = Zs;
    ZsAll(:,nn) = stripe(nn).ZRotatedSum(1:L);
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
title('local phase');

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
xlim([min(X) max(X)]);
title('unwrapped local phase');


% Image
subplot(2,4,3);
imagesc(X,1:L,ZsAll)
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
title('raw sum profile');


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
title('feedback variable');


% Wavelength
subplot(2,4,5);
plot(X,Lambda,'ko','markerfacecolor',co(1,:),'linewidth',1);
xlabel(xVar,'interpreter','none');
if isequal(xVar,'ExecutionDate')
    datetick x
end
ylabel('\lambda (sites)');
set(gca,'box','on','linewidth',1,'fontsize',10);
xlim([min(X) max(X)]);
grid on
title('wavelength');

% Rotation
subplot(2,4,6);
plot(X,Theta/pi*180,'ko','markerfacecolor',co(1,:),'linewidth',1);
xlabel(xVar,'interpreter','none');
if isequal(xVar,'ExecutionDate')
    datetick x
end
ylabel('\theta (deg)');
set(gca,'box','on','linewidth',1,'fontsize',10);
xlim([min(X) max(X)]);
grid on
title('rotation');

% Radius
subplot(2,4,7);
plot(X,Radius,'ko','markerfacecolor',co(1,:),'linewidth',1);
xlabel(xVar,'interpreter','none');
if isequal(xVar,'ExecutionDate')
    datetick x
end
ylabel('R (sites)');
set(gca,'box','on','linewidth',1,'fontsize',10);
xlim([min(X) max(X)]);
grid on
title('radius of curvature');


% Phase Histogram
subplot(2,4,8);
histogram(Phi/(2*pi),linspace(-1,1,25)*.5);
xlabel('\phi (2\pi)');
ylabel('occurences');
hold on
grid on
set(gca,'box','on','linewidth',1,'fontsize',10);
xlim([-.5 .5]);
title('phase histogram');

end

