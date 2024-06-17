function [hF,out] = dig_radialAnalysis(digdata,opts)

    if nargin == 1
        opts = struct;
    end
    
    if ~isfield(opts,'RemoveOutliers')
        opts.RemoveOutliers = 1;
    end
    
    
    %% Remove outliers
    Z   = [digdata.Zdig];
    P   = [digdata.Params];
    Xc  = [digdata.Xc_site];
    Yc  = [digdata.Yc_site];
    
    
    if opts.RemoveOutliers
        [Natoms,bad_inds] = rmoutliers([digdata.Natoms]);    
        Z(:,:,bad_inds) = [];
        P(bad_inds) = [];
        Xc(bad_inds)=[];
        Yc(bad_inds)=[];
    end
    %% Find Center
    
    
    Xcbar = round(mean(Xc));
    Ycbar = round(mean(Yc));
    
    n1 = digdata.n1;
    n2 = digdata.n2;    
    
    nlimits = [min(n1) max(n1) min(n2) max(n2)];
    %% Compute radial average
    
    L = min(abs(nlimits - [Xcbar Xcbar Ycbar Ycbar]));
    L = L-2;
    
    r = [Xcbar Xcbar Ycbar Ycbar]+[-1 1 -1 1]*L;
    
    % Indeces of bounds
    ii = [find(n1 == r(1),1) find(n1 == r(2),1) find(n2 == r(3),1) find(n2 == r(4),1)];
    
    
    % Zrad 
    Npics = size(Z,3);
    
    Zsub = zeros(length(ii(3):ii(4)),length(ii(1):ii(2)));
    for nn = 1:Npics
        Zsub(:,:,nn) = Z(ii(3):ii(4),ii(1):ii(2),nn);
    end
    
    
    % Look at average
    ZsubBar = mean(Zsub,3);
    [Tics,Average,dev,n]=radial_profile(ZsubBar,3);
    
    % Look at deviations
    ZsubDev = std(Zsub,0,3);
    [TicsDev,AverageDev,devDev,nDev]=radial_profile(ZsubDev,3);
    
    
    %% Create radial potential vector
    
    % Constants
    h = 6.62607015*10^-34;
    m = 39.96399848*1.66053906660*10^-27;
    lam = (1054*10^-9*1054*10^-9*1064*10^-9)^(1/3);
    aL = lam/2;
    
    % 120mW XDT + 2.5ER request lattice trap frequencies - calibrated
    % 02/29/24 (reanalyzed 04/11/24)
    omega_x = 2*pi*67.3;
    omega_y = 2*pi*60.2;
    omega_bar = (omega_x*omega_y).^(1/2);
    
    % Harmonic potential in Hz
    PotentialVector = 0.5*m*omega_bar^2*aL^2.*Tics.^2/h;
    
    %% Assign to output
    
    out = struct;
    out.CroppedImages                  = Zsub;
    out.AverageImage                   = ZsubBar;
    out.SiteVector1                    = r(1):r(2);
    out.SiteVector2                    = r(3):r(4);
    out.RadialCenter                   = [Xcbar Ycbar];
    out.RadialVector                   = Tics;
    out.PotentialVector                = PotentialVector;
    out.AverageOccupation              = Average;
    out.AverageOccupationUncertainty   = dev./sqrt(n);
    out.DeviationOccupation            = AverageDev;
    out.DeviationOccupationUncertainty = devDev./sqrt(nDev);


%% Plotting

    hF = figure;
    hF.Color='w';
    hF.Position = [300 100 700 500];
    clf

    if isfield(opts,'FigLabel') && ~isempty(opts.FigLabel)
        tFig=uicontrol('style','text','string',opts.FigLabel,...
            'units','pixels','backgroundcolor',...
            'w','horizontalalignment','left');
        tFig.Position(4)=tFig.Extent(4);
        tFig.Position(3)=hF.Position(3);
        tFig.Position(1:2)=[5 hF.Position(4)-tFig.Position(4)];
    end    

    % Plot radial average ndet
    subplot(221);
    errorbar(Tics(2:end),Average(2:end),dev(2:end)./sqrt(n(2:end)),'ko','markerfacecolor',[.5 .5 .5],...
        'markersize',10,'linewidth',1);
    xlabel('radial distance (sites)');
    ylabel('mean $N(r)$','interpreter','latex','fontsize',14);
    
    % Plot 2D ndet
    subplot(222);
    imagesc(r(1):r(2),r(3):r(4),ZsubBar)
    xlabel('x (sites)');
    ylabel('y (sites)');
    str = [num2str(Npics) ' images'];
    text(1,1,str,'units','pixels','fontsize',8,...
        'horizontalalignment','left','verticalalignment','bottom',...
        'color','r');
    colormap bone
    axis equal tight
    cc1=colorbar;
    cc1.Label.String = 'average of occupation';

    % Plot radial average ndet std
    subplot(223);
    plot(TicsDev(2:end),AverageDev(2:end),'ko','markerfacecolor',[.5 .5 .5],...
        'markersize',10,'linewidth',1);
    hold on;
    plot(Tics(2:end),sqrt(Average(2:end).*(1-Average(2:end))),'k.','markerfacecolor',[.25 .25 .25],...
        'markersize',10,'linewidth',1);
    xlabel('radial distance (sites)');
    ylabel('std $N(r)$','interpreter','latex','fontsize',14);
    % ylim([0 1])

    % Plot 2D ndet std average
    subplot(224);
    imagesc(r(1):r(2),r(3):r(4),ZsubDev)
    xlabel('x (sites)');
    ylabel('y (sites)');
    colormap bone
    axis equal tight
    cc=colorbar;
    cc.Label.String = 'deviation of occupation';
    % caxis(ax2.CLim);
    text(1,1,str,'units','pixels','fontsize',8,...
        'horizontalalignment','left','verticalalignment','bottom',...
        'color','r');
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
