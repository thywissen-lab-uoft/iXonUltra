function hF = dig_showFidelityMap(digdata,opts)

if nargin~=2
    opts=struct;
end

if ~isfield(opts,'doAverage')
    opts.doAverage = 1;
end

n1 = digdata.n1;
n2 = digdata.n2;

% Average of each image
z1bar = mean(digdata.Zdig(:,:,:,1),3);
n1bar = mean(sum(digdata.Zdig(:,:,:,1),[1 2]),3);
n1sigma = std(sum(digdata.Zdig(:,:,:,1),[1 2]));
z2bar = mean(digdata.Zdig(:,:,:,2),3);
n2bar = mean(sum(digdata.Zdig(:,:,:,2),[1 2]),3);
n2sigma = std(sum(digdata.Zdig(:,:,:,1),[1 2]));

lostbarN = mean(digdata.lost_number);
lostbarP = mean(digdata.lost_fraction);
hopbarN = mean(digdata.hop_number);
hopbarP = mean(digdata.hop_fraction);
% Make the figure
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
colormap(ax1,'bone')
set(ax1,'fontsize',10,'fontname','times','box','on',...
    'linewidth',1,'ydir','normal');    
title_str = ['image 1 : $\langle N \rangle=' num2str(n1bar) '$'];
title(ax1,title_str,'interpreter','latex');    
axis(ax1,'tight');
axis(ax1,'equal');
cL=get(ax1,'CLim');
colorbar(ax1);
xlabel(ax1,'site 1')
ylabel(ax1,'site 2');

% Image 2
ax2=subplot(2,2,2,'parent',hF);    
imagesc(n1,n2,z2bar,'parent',ax2);
colormap(ax2,'bone')
set(ax2,'fontsize',10,'fontname','times','box','on',...
    'linewidth',1,'ydir','normal');    
title_str = ['image 2 : $\langle N \rangle=' num2str(n2bar) '$'];
title(ax2,title_str,'interpreter','latex');    
axis(ax2,'tight');
axis(ax2,'equal');
colorbar(ax2);
xlabel(ax2,'site 1')
ylabel(ax2,'site 2');
caxis(ax2,cL);

% Differential Image
ax3=subplot(2,2,3,'parent',hF);    
imagesc(n1,n2,z1bar-z2bar,'parent',ax3);
colormap(ax3,'bone')
s3 = ['$\mathrm{lost}:' num2str(lostbarN) ...
    '~(' num2str(round(lostbarP*100,1)) ' \%),~'  ...
    '\mathrm{hop}:' num2str(hopbarN) ...
    '~(' num2str(round(hopbarP*100,1)) ' \%)$'];
% caxis([-1 1]);
colorbar(ax3)
xlabel(ax3,'site 1')
ylabel(ax3,'site 2');
set(ax3,'fontsize',10,'fontname','times','box','on','linewidth',1,'ydir','normal');
title(ax3,'image 1 - image 2');  
axis(ax3,'tight');
axis(ax3,'equal');
text(.01,.01,s3,'units','normalized','fontsize',8,'color','black',...
    'verticalalignment','bottom','horizontalalignment','left',...
    'interpreter','latex','backgroundcolor',[1 1 1 .8],'margin',1,...
    'parent',ax3)

    end



