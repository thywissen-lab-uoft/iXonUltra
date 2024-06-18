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
    
    %% Recenter every image to have the same mean?
    
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

    % Average Radial Deviation profile (and std of std)
    radial_charge_std_mean  = mean(radial_charge_err,2);
    radial_charge_std_std  = std(radial_charge_err,1,2);
    
    
    

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
    out.AverageOccupationUncertainty   = radial_charge_mean_dev./sqrt(n);
    out.DeviationOccupation            = radial_charge_std_mean;
%     out.DeviationOccupationUncertainty = devDev./sqrt(nDev);

%% Plotting

    hF = figure;
    hF.Color='w';
    hF.Position = [300 100 1500 350];
    clf

    if isfield(opts,'FigLabel') && ~isempty(opts.FigLabel)
        tFig=uicontrol('style','text','string',opts.FigLabel,...
            'units','pixels','backgroundcolor',...
            'w','horizontalalignment','left');
        tFig.Position(4)=tFig.Extent(4);
        tFig.Position(3)=300;
        tFig.Position(1:2)=[5 hF.Position(4)-tFig.Position(4)];
    end    
    
    % 2D Image of Average Image
    subplot(131);
    imagesc(r(1):r(2),r(3):r(4),imboxfilt(ZsubBar,opts.BinStep))
    xlabel('x (sites)');
    ylabel('y (sites)');
    
    strSummary = ['$n=' num2str(Npics) ' $ images' newline ...
        '$N=' num2str(round(mean(digdata.Natoms))) '\pm' num2str(round(std(digdata.Natoms,1))) '$ atoms'];
    text(1,1,strSummary,'units','pixels','fontsize',12,...
        'horizontalalignment','left','verticalalignment','bottom',...
        'color','r','interpreter','latex');  
    
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
    subplot(132);
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

    % Plot radial average ndet std
    subplot(133);        
    plot(rVec(2:end),radial_charge_mean_std(2:end)/sqrt(Npics),'ko','markerfacecolor',[.5 .5 .5],...
        'markersize',8,'linewidth',1);
    title('uncertainty of average');
    ylabel('std $\{N(r)\}/\sqrt{n}$','interpreter','latex','fontsize',14);
    xlabel('radial distance (sites)');    
    text(.99,.99,strSummary,'units','normalized','fontsize',12,'horizontalalignment','right',...
        'verticalalignment','top','interpreter','latex');
    
    if isfield(opts,'rMaxShow')
       xlim([0 opts.rMaxShow]); 
    end


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
