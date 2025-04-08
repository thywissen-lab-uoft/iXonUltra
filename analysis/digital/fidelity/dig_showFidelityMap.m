function hF = dig_showFidelityMap(digdata,opts)

if nargin~=2
    opts=struct;
end

if ~isfield(opts,'doAverage')
    opts.doAverage = 1;
end


%% Calculate Stuff
n1 = digdata.n1;
n2 = digdata.n2;

% Average of each image
z1bar       = mean(digdata.Zdig(:,:,:,1),3);
n1bar       = mean(sum(digdata.Zdig(:,:,:,1),[1 2]),3);
n1sigma     = std(sum(digdata.Zdig(:,:,:,1),[1 2]));
z2bar       = mean(digdata.Zdig(:,:,:,2),3);
n2bar       = mean(sum(digdata.Zdig(:,:,:,2),[1 2]),3);
n2sigma     = std(sum(digdata.Zdig(:,:,:,2),[1 2]));

lostbarN    = mean(digdata.lost_number);
lostbarP    = mean(digdata.lost_fraction);
hopbarN     = mean(digdata.hop_number);
hopbarP     = mean(digdata.hop_fraction);


%% Calculte Histgoram

[nn1, nn2]=meshgrid(n1,n2);  
dzbar = z1bar-z2bar;
zbar = 0.5*(z1bar+z2bar);

n1c = sum(zbar.*nn1,'all')/sum(zbar,'all');
n2c = sum(zbar.*nn2,'all')/sum(zbar,'all');

n1_variance = sum(zbar.*(nn1-n1c).^2,'all')/sum(zbar,'all');
n2_variance = sum(zbar.*(nn1-n2c).^2,'all')/sum(zbar,'all');

s1 = sqrt(n1_variance);
s2 = sqrt(n2_variance);

sigma = sqrt(s2*s1);

nnr = sqrt((nn1-n1c).^2+(nn2-n2c).^2);

r_hist_defect = [];
r_hist_atom = [];
for kk=1:size(digdata.Zdig,3)
    z1 = logical(digdata.Zdig(:,:,kk,1));
    dz = logical(abs(digdata.Zdig(:,:,kk,1)-digdata.Zdig(:,:,kk,2)));
    r_this_defect = nnr(dz);
    r_this_defect = r_this_defect(:);
    r_this_atom = nnr(z1);
    r_this_atom = r_this_atom(:);
    r_hist_defect=[r_hist_defect; r_this_defect];
    r_hist_atom=[r_hist_atom; r_this_atom];
end
edges=linspace(0,60,20);   

[N_defect,edges_defect]=histcounts(r_hist_defect(:),edges);
[N_atom,edges_atom]=histcounts(r_hist_atom(:),edges);
centers = (edges(1:end-1) + edges(2:end))/2;   


bad_inds = [N_atom<5]

%% Make the figure
if ~isfield(opts,'Parent') || isempty(opts.Parent)
    opts.Parent = figure;
    set(opts.Parent,'color','w','Name',['fidelity_map']);
    clf
    opts.Parent.Position=[0 50 900 700];
end

hF=opts.Parent;

if isfield(opts,'FigLabel') && ~isempty(opts.FigLabel)
    t=uicontrol('style','text','string',opts.FigLabel,'fontsize',7,...
        'backgroundcolor','w','Position',[1 1 600 15],'horizontalalignment','left',...
        'parent',hF);
    t.Position(2) = t.Parent.Position(4)-t.Position(4)-2;
end
    
% Image 1
ax1=subplot(2,2,1,'parent',hF);    
imagesc(n1,n2,z1bar,'parent',ax1);
colormap(ax1,slanCM('magma'))
set(ax1,'fontsize',10,'fontname','times','box','on',...
    'linewidth',1,'ydir','normal');    
axis(ax1,'tight');
axis(ax1,'equal');
cL=get(ax1,'CLim');
caxis(ax1,[0 cL(2)]);
colorbar(ax1);
xlabel(ax1,'site 1')
ylabel(ax1,'site 2');
if size(digdata.Zdig,4)>1
    s1 = ['$\langle N \rangle = ' num2str(round(n1bar)) ...
        '\pm' num2str(round(n1sigma)) '$'];
        tit_str=['mean(image 1)'];
else
    s1 = ['$\langle N \rangle = '  num2str(round(n1bar)) '$'];
    tit_str=['image 2'];
end
text(.01,.01,s1,'units','normalized','fontsize',8,'color','black',...
    'verticalalignment','bottom','horizontalalignment','left',...
    'interpreter','latex','backgroundcolor',[1 1 1 .8],'margin',1,...
    'parent',ax1)
title_str = tit_str;
title(ax1,title_str,'interpreter','latex');    
xlim(ax1,n1c+[-1 1]*sigma*2.5);
ylim(ax1,n2c+[-1 1]*sigma*2.5);


% Image 2
ax2=subplot(2,2,2,'parent',hF);    
imagesc(n1,n2,z2bar,'parent',ax2);
colormap(ax2,slanCM('magma'))
set(ax2,'fontsize',10,'fontname','times','box','on',...
    'linewidth',1,'ydir','normal');    
axis(ax2,'tight');
axis(ax2,'equal');
colorbar(ax2);
xlabel(ax2,'site 1')
ylabel(ax2,'site 2');
caxis(ax2,[0 cL(2)]);
if size(digdata.Zdig,4)>1
    s2 = ['$\langle N \rangle = ' num2str(round(n2bar)) ...
        '\pm' num2str(round(n2sigma)) '$'];
    tit_str=['mean(image 2)'];
else
    s2 = ['$\langle N \rangle = '  num2str(round(n2bar)) '$'];
    tit_str=['image 2'];
end
text(.01,.01,s2,'units','normalized','fontsize',8,'color','black',...
    'verticalalignment','bottom','horizontalalignment','left',...
    'interpreter','latex','backgroundcolor',[1 1 1 .8],'margin',1,...
    'parent',ax2)
title_str = tit_str;
title(ax2,title_str,'interpreter','latex');    
xlim(ax2,n1c+[-1 1]*sigma*2.5);
ylim(ax2,n2c+[-1 1]*sigma*2.5);

% Differential Image
ax3=subplot(2,2,3,'parent',hF);    
cla(ax3);
imagesc(n1,n2,z1bar-z2bar,'parent',ax3);
% colormap(ax3,'bone')
colormap(ax3,blue_black_red)
s3 = ['$\mathrm{lost}:' num2str(round(lostbarN)) ...
    '~(' num2str(round(lostbarP*100,1)) ' \%),~'  ...
    '\mathrm{hop}:' num2str(round(hopbarN)) ...
    '~(' num2str(round(hopbarP*100,1)) ' \%)$'];
colorbar(ax3)
xlabel(ax3,'site 1')
ylabel(ax3,'site 2');
set(ax3,'fontsize',10,'fontname','times','box','on','linewidth',1,'ydir','normal');
axis(ax3,'tight');
axis(ax3,'equal');
text(.01,.01,s3,'units','normalized','fontsize',8,'color','black',...
    'verticalalignment','bottom','horizontalalignment','left',...
    'interpreter','latex','backgroundcolor',[1 1 1 .8],'margin',1,...
    'parent',ax3)
cLD = max(abs(get(ax3,'CLim')));
caxis(ax3,[-1 1]*cLD)

if size(digdata.Zdig,4)>1    
    title_str=['mean(image 1 - image 2)'];    
else 
    title_str=['image 1 - image 2'];
end
title(ax3,title_str,'interpreter','latex');  
xlim(ax3,n1c+[-1 1]*sigma*2.5);
ylim(ax3,n2c+[-1 1]*sigma*2.5);


    ax4 = subplot(4,4,[11 12]);
    cla(ax4);
    yyaxis left
    bar(centers,N_atom,'parent',ax4);
    ylabel(ax4,'N atom(r)')
    drawnow;
    title(ax4,'atom histogram','interpreter','latex')

    ax5 = subplot(4,4,[15 16]);
    cla(ax5);
    yyaxis left
    bar(centers,N_defect,'parent',ax5);
    ylabel(ax5,'N defect(r)')
    ylim(ax5,get(ax4,'YLim'))
    yyaxis right
    plot(centers(~bad_inds),N_defect(~bad_inds)./N_atom(~bad_inds),'parent',ax5);
    ylabel(ax5,'defect per atom (r)')
    title(ax5,'defect histogram','interpreter','latex')
    xlabel(ax5,'radius (sites)')

end

function cc=blue_black_red(N)
if nargin==0
N=256;
end

black2red=[linspace(0,1,N);linspace(0,0,N);linspace(0,0,N)]';
blue2black=[linspace(0,0,N);linspace(0,0,N);linspace(1,0,N)]';

cc = [blue2black;black2red];

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

