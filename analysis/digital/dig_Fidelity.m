function out = dig_Fidelity(Zdig,n1,n2,opts)

if nargin <4
    opts=struct;
    opts.FigureNumber = 4001;
end
img1 = Zdig(:,:,1);
img2 = Zdig(:,:,2);
dImg = img2-img1;

N1 = sum(img1,'all');
N2 = sum(img2,'all');

% Number of lost atoms
Nlost = N1-N2;
Nlost_percent = Nlost/N1;

% Number of atoms that have hopped
Nhop = sum((img2-img1)>0,'all'); % look for "new" atoms in image 2
Nhop_percent = Nhop/N1;

%% Radial

Zmap = (img1+img2);
Zmap = Zmap/sum(Zmap,'all');

[xx,yy]=meshgrid(n1,n2);

xc = round(sum(xx.*Zmap,'all'));
yc = round(sum(yy.*Zmap,'all'));

R = sqrt((xx-xc).^2+(yy-yc).^2); % distance from center
iEvent = logical(abs(dImg));

Revent = R(iEvent);


nlimits = [min(n1) max(n1) min(n2) max(n2)];    
L = min(abs(nlimits - [xc xc yc yc]));
L = L-2;
r = [xc xc yc yc]+[-1 1 -1 1]*L;

% Indeces of bounds
ii = [find(n1 == r(1),1) find(n1 == r(2),1) find(n2 == r(3),1) find(n2 == r(4),1)];
img_avg = 0.5*(img1+img2);
img_avg_sub = img_avg(ii(3):ii(4),ii(1):ii(2));
dR = 10;

[r,N_expect,dev,n]=radial_profile(img_avg_sub,dR);

edges=0:dR:max(r);

%%
hF = figure(opts.FigureNumber) ;
set(hF,'color','w','Name','Dig Fidelity');
clf
hF.Position=[0 710 1100 250];

if isfield(opts,'Label') && ~isempty(opts.Label)
uicontrol('style','text','string',opts.Label,'fontsize',7,...
    'backgroundcolor','w','Position',[1 1 600 15],'horizontalalignment','center');
end
% subplot(221);
subplot(141);

imagesc(n1,n2,img1);
colormap bone
caxis([0 1]);
xlabel('site 1')
ylabel('site 2');
set(gca,'fontsize',10,'fontname','times','box','on','linewidth',1,'ydir','normal');

title_str = ['image 1 : $N=' num2str(N1) '$'];
title(title_str,'interpreter','latex');

axis equal tight

hold on
plot(xc,yc,'o','color','r','markersize',5,'markerfacecolor','r');


% subplot(222);
subplot(142);

imagesc(n1,n2,img2);
colormap bone
caxis([0 1]);
xlabel('site 1')
ylabel('site 2');
set(gca,'fontsize',10,'fontname','times','box','on','linewidth',1,'ydir','normal');
title_str = ['image 1 : $N=' num2str(N2) '$'];
title(title_str,'interpreter','latex');
axis equal tight

hold on
plot(xc,yc,'o','color','r','markersize',5,'markerfacecolor','r');

% subplot(223);
subplot(143);

co=get(gca,'colororder');
imagesc(n1,n2,dImg);
colormap bone
s3 = ['$\mathrm{lost}:' num2str(Nlost) '~(' num2str(round(Nlost_percent*100,1)) ' \%),~'  ...
    '\mathrm{hop}:' num2str(Nhop) '~(' num2str(round(Nhop_percent*100,1)) ' \%)$'];
caxis([-1 1]);
xlabel('site 1')
ylabel('site 2');
set(gca,'fontsize',10,'fontname','times','box','on','linewidth',1,'ydir','normal');
title('image 2 - image1');



axis equal tight
text(.01,.01,s3,'units','normalized','fontsize',8,'color','black',...
    'verticalalignment','bottom','horizontalalignment','left',...
    'interpreter','latex','backgroundcolor',[1 1 1 .8],'margin',1)
hold on
plot(xc,yc,'o','color','r','markersize',3,'markerfacecolor','r');


cc = parula(length(edges));
N = length(edges);

% rList = [20 40 60 80];
tt=linspace(0,2*pi,100);
for kk=2:length(edges)
    plot(xc+edges(kk)*cos(tt),yc+edges(kk)*sin(tt),'color',cc(kk,:))
end
xlim([min(n1) max(n1)]);
ylim([min(n2) max(n2)]);

% subplot(4,4,[11 12]);
subplot(2,4,[4]);

histogram(Revent,edges);
xlabel('radial position (sites)')
ylabel('loss or hop occurence')
xlim([0 max(edges)]);
set(gca,'box','on','linewidth',1,'fontname','times','fontsize',8);

% subplot(4,4,[15 16]);
subplot(2,4,[8]);

co=get(gca,'colororder');
xlabel('radial position (sites)')
% plot(r,N_expect,dev,'o','markerfacecolor',co(1,:),'linewidth',1,'markersize',6,...
    % 'color',co(1,:)*.5);
    % keyboard
    plot(r(1:N),N_expect(1:N),'k-');
    hold on

ps=scatter(r(1:N),N_expect(1:N),40,cc,'filled');
set(ps,'markeredgecolor','k')
ylabel('average occupation')
xlim([0 max(edges)]);
xlabel('radial position (sites)')
set(gca,'box','on','linewidth',1,'fontname','times','fontsize',8);
%%
out = struct;

end

function [Tics,Average,dev,n]=radial_profile(data,radial_step)
%main axii cpecified:
x=(1:size(data,2))-size(data,2)/2;
y=(1:size(data,1))-size(data,1)/2;
% coordinate grid:
[X,Y]=meshgrid(x,y);
% creating circular layers
Z_integer=round(abs(X+1i*Y)/radial_step)+1;
% % illustrating the principle:
% % figure;imagesc(Z_integer.*data)
% very fast MatLab calculations:
Tics=accumarray(Z_integer(:),abs(X(:)+1i*Y(:)),[],@mean);
Average=accumarray(Z_integer(:),data(:),[],@mean);


dev=accumarray(Z_integer(:),data(:),[],@std);

n= accumarray(Z_integer(:),data(:),[],@(x) numel(x));

end