% function out = dig_Fidelity(Zdig,n1,n2,opts)
function [fidelity,hF_out,hF] = dig_Fidelity(digdata,opts)
% Author : CJ Fujiwara
%
% This codes calculates the fidelity of digitized data
hF_out=[];
hF=[];

%% Settings
if nargin ==1
    opts = struct;
end

if length(digdata.FileNames)==1
    opts.doDebug = 1;
else
    opts.doDebug = 0;
end

opts.FigureNumber = 4001;
%% Grab data
  
n1 = digdata.n1;
n2 = digdata.n2;

% z1_all= zeros(length(n1),length(n2),length(digdata.FileNames));
fidelity=struct;

%% Iterate
for kk = 1:length(digdata.FileNames)
    z1 = digdata.Zdig(:,:,kk,1);    % Image 1
    z2 = digdata.Zdig(:,:,kk,2);    % Image 2

    N1 = sum(z1,'all');           % Atom Number 1
    N2 = sum(z2,'all');           % Atom Number 2

    dz = z1 - z2;                   % Differential Image

    Nlost = N1 - N2;                % Number Lost
    Rlost = Nlost/N1;               % Rate of Loss

    Nhop = sum(abs(dz),'all')-Nlost;% Number Hopped
    Rhop = Nhop/N1;              % Rate of hop
  
    %% Radial
    
    % Zmap = (img1+img2);
    % Zmap = Zmap/sum(Zmap,'all');
    % 
    % [xx,yy]=meshgrid(n1,n2);
    % 
    % % keyboard
    % 
    % xc = round(sum(xx.*Zmap,'all'));
    % yc = round(sum(yy.*Zmap,'all'));
    % 
    % R = sqrt((xx-xc).^2+(yy-yc).^2); % distance from center
    % iEvent = logical(abs(dImg));
    % 
    % Revent = R(iEvent);
    % 
    % 
    % nlimits = [min(n1) max(n1) min(n2) max(n2)];    
    % L = min(abs(nlimits - [xc xc yc yc]));
    % L = L-2;
    % r = [xc xc yc yc]+[-1 1 -1 1]*L;
    % 
    % % Indeces of bounds
    % ii = [find(n1 == r(1),1) find(n1 == r(2),1) find(n2 == r(3),1) find(n2 == r(4),1)];
    % img_avg = 0.5*(img1+img2);
    % img_avg_sub = img_avg(ii(3):ii(4),ii(1):ii(2));
    % dR = 10;
    % 
    % [r,N_expect,dev,n]=radial_profile(img_avg_sub,dR);
    % 
    % edges=0:dR:max(r);
    
    %%
    if opts.doDebug
        hF = figure(opts.FigureNumber+kk) ;
        set(hF,'color','w','Name',['fidelity_' num2str(kk)]);
        clf
        hF.Position=[0 710 1100 250];
    if isfield(opts,'FigLabel') && ~isempty(opts.FigLabel)
        t=uicontrol('style','text','string',opts.FigLabel,'fontsize',7,...
            'backgroundcolor','w','Position',[1 1 600 15],'horizontalalignment','left');
        t.Position(2) = t.Parent.Position(4)-t.Position(4)-2;

    end
    
        % Image 1
        subplot(131);    
        imagesc(n1,n2,z1);
        colormap jet
        caxis([0 1]);
        set(gca,'fontsize',10,'fontname','times','box','on',...
            'linewidth',1,'ydir','normal');    
        title_str = ['image 1 : $N=' num2str(N1) '$'];
        title(title_str,'interpreter','latex');    
        axis equal tight    
        
        % Image 2
        subplot(132);    
        imagesc(n1,n2,z2);
        colormap jet
        caxis([0 1]);
        set(gca,'fontsize',10,'fontname','times','box','on',...
            'linewidth',1,'ydir','normal');
        title_str = ['image 1 : $N=' num2str(N2) '$'];
        title(title_str,'interpreter','latex');
        axis equal tight    
        
        % Differential Image
        subplot(133);    
        imagesc(n1,n2,dz);
        colormap jet
        s3 = ['$\mathrm{lost}:' num2str(Nlost) ...
            '~(' num2str(round(Rlost*100,1)) ' \%),~'  ...
            '\mathrm{hop}:' num2str(Nhop) ...
            '~(' num2str(round(Rhop*100,1)) ' \%)$'];
        caxis([-1 1]);
        xlabel('site 1')
        ylabel('site 2');
        set(gca,'fontsize',10,'fontname','times','box','on','linewidth',1,'ydir','normal');
        title('image 1 - image 2');  
        axis equal tight
        text(.01,.01,s3,'units','normalized','fontsize',8,'color','black',...
            'verticalalignment','bottom','horizontalalignment','left',...
            'interpreter','latex','backgroundcolor',[1 1 1 .8],'margin',1)
        hold on
    end
    % plot(xc,yc,'o','color','r','markersize',3,'markerfacecolor','r');
    
    % 
    % cc = parula(length(edges));
    % N = length(edges);
    % 
    % % rList = [20 40 60 80];
    % tt=linspace(0,2*pi,100);
    % for kk=2:length(edges)
    %     plot(xc+edges(kk)*cos(tt),yc+edges(kk)*sin(tt),'color',cc(kk,:))
    % end
    % xlim([min(n1) max(n1)]);
    % ylim([min(n2) max(n2)]);
    % 
    % % subplot(4,4,[11 12]);
    % subplot(2,4,[4]);
    % 
    % histogram(Revent,edges);
    % xlabel('radial position (sites)')
    % ylabel('loss or hop occurence')
    % xlim([0 max(edges)]);
    % set(gca,'box','on','linewidth',1,'fontname','times','fontsize',8);
    % 
    % % subplot(4,4,[15 16]);
    % subplot(2,4,[8]);
    % 
    % co=get(gca,'colororder');
    % xlabel('radial position (sites)')
    % % plot(r,N_expect,dev,'o','markerfacecolor',co(1,:),'linewidth',1,'markersize',6,...
    %     % 'color',co(1,:)*.5);
    %     % keyboard
    %     plot(r(1:N),N_expect(1:N),'k-');
    %     hold on
    % 
    % ps=scatter(r(1:N),N_expect(1:N),40,cc,'filled');
    % set(ps,'markeredgecolor','k')
    % ylabel('average occupation')
    % xlim([0 max(edges)]);
    % xlabel('radial position (sites)')
    % set(gca,'box','on','linewidth',1,'fontname','times','fontsize',8);
    %%
    fidelity(kk).N1          = N1;
    fidelity(kk).N2          = N2;
    fidelity(kk).Nlost       = Nlost;
    fidelity(kk).Rlost       = Rlost;
    fidelity(kk).Nhop        = Nhop;
    fidelity(kk).Rhop        = Rhop;
end

digdata.Fidelity = fidelity;

%% Summary Figure

if length(digdata.FileNames)>1
    hF_out = figure;
    set(hF_out,'color','w','Name',['fidelity_summary']);
    clf
    hF_out.Position=[0 710 1300 250];
    if isfield(opts,'FigLabel') && ~isempty(opts.FigLabel)
        t=uicontrol('style','text','string',opts.FigLabel,'fontsize',7,...
            'backgroundcolor','w','Position',[1 1 600 15],'horizontalalignment','left');
        t.Position(2) = t.Parent.Position(4)-t.Position(4)-2;
    end

    X= [digdata.X];
    co=get(gca,'colororder');
    subplot(141);
    p1=plot(X,[fidelity.N1],'ko','markerfacecolor',co(4,:),'color',...
        co(4,:)*.5,'linewidth',1);
    hold on
    p2=plot(X,[fidelity.N2],'^','markerfacecolor',co(5,:),'color',...
        co(5,:)*.5,'linewidth',1);
    ylabel('atom number')
    if isequal(digdata.xVar,'ExecutionDate')
        datetick x
    end
    legend([p1 p2],{'N1','N2'});
    yL=get(gca,'YLim');
    set(gca,'YLim',[0 yL(2)]);
    xlabel(digdata.xVar,'interpreter','none');

    subplot(142);
    pLost=plot(X,[fidelity.Nlost],'o','markerfacecolor',co(2,:),'color',...
        co(2,:)*.5,'linewidth',1);
    hold on
    pHop=plot(X,[fidelity.Nhop],'o','markerfacecolor',co(1,:),'color',...
        co(1,:)*.5,'linewidth',1);
    legend([pLost pHop],{'N lost','N hop'});
    ylabel('number')
    if isequal(digdata.xVar,'ExecutionDate')
        datetick x
    end
    yL=get(gca,'YLim');
    set(gca,'YLim',[0 yL(2)]);
    xlabel(digdata.xVar,'interpreter','none');

    subplot(143);
    pRLost=plot(X,[fidelity.Rlost],'o','markerfacecolor',co(2,:),'color',...
        co(2,:)*.5,'linewidth',1);
    hold on
    pRHop=plot(X,[fidelity.Rhop],'o','markerfacecolor',co(1,:),'color',...
        co(1,:)*.5,'linewidth',1);
    legend([pRLost pRHop],{'rate lost','rate hop'});
    ylabel('number')
    if isequal(digdata.xVar,'ExecutionDate')
        datetick x
    end
    yL=get(gca,'YLim');
    set(gca,'YLim',[0 yL(2)]);
    xlabel(digdata.xVar,'interpreter','none');

    subplot(144);
    dz = digdata.Zdig(:,:,:,1)-digdata.Zdig(:,:,:,2);
    defects = abs(dz);
    defects = sum(defects,3);
    imagesc(n1,n2,defects);
    colormap jet
    title('defect map');
    colorbar 
    axis equal tight

    xc = mean(digdata.Xc_site,'all');
    yc = mean(digdata.Yc_site,'all');

    xlim(xc+[-50 50])
    ylim(yc+[-50 50])


end

end

% function [Tics,Average,dev,n]=radial_profile(data,radial_step)
% %main axii cpecified:
% x=(1:size(data,2))-size(data,2)/2;
% y=(1:size(data,1))-size(data,1)/2;
% % coordinate grid:
% [X,Y]=meshgrid(x,y);
% % creating circular layers
% Z_integer=round(abs(X+1i*Y)/radial_step)+1;
% % % illustrating the principle:
% % % figure;imagesc(Z_integer.*data)
% % very fast MatLab calculations:
% Tics=accumarray(Z_integer(:),abs(X(:)+1i*Y(:)),[],@mean);
% Average=accumarray(Z_integer(:),data(:),[],@mean);
% 
% 
% dev=accumarray(Z_integer(:),data(:),[],@std);
% 
% n= accumarray(Z_integer(:),data(:),[],@(x) numel(x));
% 
% end