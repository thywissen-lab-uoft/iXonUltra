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
    
    %% Compute Square ROI around center of cloud for radial computing    
    L = min(abs(nlimits - [Xcbar Xcbar Ycbar Ycbar]));
    L = L-2;    
    r = [Xcbar Xcbar Ycbar Ycbar]+[-1 1 -1 1]*L;
    
    % Indeces of bounds
    ii = [find(n1 == r(1),1) find(n1 == r(2),1) find(n2 == r(3),1) find(n2 == r(4),1)];
    
    %% Compute Average Image
    
    Npics = size(Z,3);    
    Zsub = zeros(length(ii(3):ii(4)),length(ii(1):ii(2)));
    
    for nn = 1:Npics
        Zsub(:,:,nn) = Z(ii(3):ii(4),ii(1):ii(2),nn);
    end    
    ZsubBar = mean(Zsub,3);
    
    % Radial profile of average image
    [rVec,radial_charge_mean,radial_charge_mean_dev,n]=radial_profile(ZsubBar,opts.BinStep);

    %% Compute Radial Profile of each Image
    radial_charge = zeros(length(rVec),Npics);      %  Radial Charge Matrix
    radial_charge_err = zeros(length(rVec),Npics);  %  Radial Charge  Deviation Matrix
    
    % Iterate over each each
    for ii = 1:Npics
        % Compute radial profile
        [rVec,charge,charge_dev,n]=...
            radial_profile(Zsub(:,:,ii),opts.BinStep);
        radial_charge(:,ii) = charge;
        radial_charge_err(:,ii) = charge_dev./sqrt(n);
    end  
    
    % Average Radial charge profile
    radial_charge_mean      = mean(radial_charge,2);     
    radial_charge_mean_std  = std(radial_charge,1,2);     
    

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
    PotentialVector = 0.5*m*omega_bar^2*aL^2.*rVec.^2/h;
    
    %% Assign to output
    
    out = struct;
    out.CroppedImages                  = Zsub;
    out.AverageImage                   = ZsubBar;
    out.SiteVector1                    = r(1):r(2);
    out.SiteVector2                    = r(3):r(4);
    out.RadialCenter                   = [Xcbar Ycbar];
    out.RadialVector                   = rVec;
    out.PotentialVector                = PotentialVector;
    out.AverageOccupation              = radial_charge_mean;
    out.AverageOccupationUncertainty   = radial_charge_mean_std/sqrt(Npics);

%% Plotting

    hF = figure;
    hF.Color='w';
    hF.Position = [300 100 1200 350];
    W= hF.Position(3);
    H = hF.Position(4);
    clf

    if isfield(opts,'FigLabel') && ~isempty(opts.FigLabel)
        tFig=uicontrol('style','text','string',opts.FigLabel,...
            'units','pixels','backgroundcolor',...
            'w','horizontalalignment','left');
        tFig.Position(4)=tFig.Extent(4);
        tFig.Position(3)=tFig.Extent(3);
%         tFig.Position(1:2)=[5 hF.Position(4)-tFig.Position(4)];
        tFig.Position(1:2)=[1 1];

    end    
    
    % 2D Image of Average Image
    ax1=subplot(131);
    ax1.Units='pixels';
    [x0, y0, w, h]=getAxesPos(1,3,1,W,H);
    ax1.Position = [x0 y0 w h];
    
    
    imagesc(r(1):r(2),r(3):r(4),imboxfilt(ZsubBar,opts.BinStep))
    xlabel('position (sites)');
    ylabel('position (sites)');
    
    strBoxCar = ['boxcar avg: ' num2str(opts.BinStep) '\times' num2str(opts.BinStep)];
    
    strSummary = ['$n=' num2str(Npics) ' $ images' newline ...
        '$N=' num2str(round(mean(digdata.Natoms))) '\pm' ...        
        num2str(round(std(digdata.Natoms,1))) '$ atoms'];    
    
    text(1,1,strSummary,'units','pixels','fontsize',10,...
        'horizontalalignment','left','verticalalignment','bottom',...
        'color','r','interpreter','latex');  
    
    text(.01,.99,strBoxCar,'units','normalized','fontsize',10,...
        'verticalalignment','top','horizontalalignment','left',...
        'color','r');
        ax1.XAxisLocation='Top';
        set(ax1,'FontSize',8);

    colormap bone
    axis equal tight
    cc1=colorbar;
    cc1.Label.String = 'charge occupation';
    title('average image');
    set(gca,'ydir','normal','box','on','linewidth',1);
    
    if isfield(opts,'nMaxShow')
       caxis([0 opts.nMaxShow]); 
    end

    % Average charge density
    ax2=subplot(132);
    ax2.Units='pixels';
    [x0, y0, w, h]=getAxesPos(1,3,2,W,H);
    ax2.Position = [x0 y0 w h];
    
    errorbar(rVec(2:end),radial_charge_mean(2:end),...
        radial_charge_mean_std(2:end)/sqrt(Npics),...
        'ko','markerfacecolor',[.5 .5 .5],...
        'markersize',8,'linewidth',1);
    hold on      
    xlabel('radial distance (sites)');
    ylabel('mean $\{N(r)\}$','interpreter','latex','fontsize',14);    
    title('average radial profile');        
    text(.99,.99,strSummary,'units','normalized','fontsize',12,'horizontalalignment','right',...
        'verticalalignment','top','interpreter','latex');
    
    if isfield(opts,'nMaxShow')
       ylim([0 opts.nMaxShow]); 
    end
    
    if isfield(opts,'rMaxShow')
       xlim([0 opts.rMaxShow]); 
    end
    
    strRadialBin = ['radial bin ' char(916) 'r:' num2str(opts.BinStep)];
    text(.01,.99,strRadialBin,'units','normalized','horizontalalignment','left',...
        'verticalalignment','top','fontsize',10);

    % Plot radial average ndet std
    ax3=subplot(133);
    ax3.Units='pixels';
    [x0, y0, w, h]=getAxesPos(1,3,3,W,H);
    ax3.Position = [x0 y0 w h];
    
    plot(rVec(2:end),radial_charge_mean_std(2:end)/sqrt(Npics),...
        'ko','markerfacecolor',[.5 .5 .5],...
        'markersize',8,'linewidth',1);
    title('standard error');
    ylabel('std $\{N(r)\}/\sqrt{n}$','interpreter','latex','fontsize',14);
    xlabel('radial distance (sites)');    
    text(.99,.99,strSummary,'units','normalized','fontsize',12,'horizontalalignment','right',...
        'verticalalignment','top','interpreter','latex');
    
    if isfield(opts,'rMaxShow')
       xlim([0 opts.rMaxShow]); 
    end
    text(.01,.99,strRadialBin,'units','normalized','horizontalalignment','left',...
        'verticalalignment','top','fontsize',10);

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


% dev=accumarray(Z_integer(:),data(:),[],@std);
dev=accumarray(Z_integer(:),data(:),[],@(x) std(x,1));

n= accumarray(Z_integer(:),data(:),[],@(x) numel(x));

end


function [axX,axY,axWidth,axHeight]=getAxesPos(num_rows,num_cols,index,figW,figH)
yTop=30;
yBot=50;

xLeft=35;
xRight=20;

ySpace=35;
xSpace=70;

rowNumber = floor((index-1)/num_cols)+1;
colNumber = mod(index-1,num_cols)+1;

axHeight = (figH-yTop-yBot-ySpace*(num_rows-1))/num_rows;
axWidth = (figW-xLeft-xRight-xSpace*(num_cols-1))/num_cols;

axX = xLeft + (axWidth+xSpace)*(colNumber-1);
axY=(figH-yTop-axHeight)-(rowNumber-1)*(axHeight+ySpace);

end

