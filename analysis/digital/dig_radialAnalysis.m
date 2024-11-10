function [hF,out] = dig_radialAnalysis(digdata,opts)
    if nargin == 1
        opts = struct;
    end    
    if ~isfield(opts,'RemoveOutliers')
        opts.RemoveOutliers = 1;
    end     

    if ~isfield(opts,'BinStep')
        opts.BinStep = 3;
    end   

    if ~isfield(opts,'TrapOmega')
        opts.TrapOmega = 2*pi*65;%omega
    end

    if ~isfield(opts,'Tunneling')
        opts.Tunneling = 563;%hz
    end

    % Constants
    h = 6.62607015e-34;         % planck's constant [Js]
    amu = 1.66053906660e-27;    % atomic mass unit[kg]    
    m = amu*39.96399848;            % 40K mass [kg]
    aL = (1054e-9)/2;           % lattice spacing [m]
    kB= 1.380649e-23;           % J/K

    site2Vr_Hz = @(n) 0.5*m*opts.TrapOmega^2*aL^2*n.^2/h;
    sigma2T_nK = @(sigmaR) 1e9*sigmaR^2*m*opts.TrapOmega^2/kB;

    nstar = linspace(0,100,1e4);
    istar=find([site2Vr_Hz(nstar)]>=4*opts.Tunneling,1);
    rstar=nstar(istar)*0.527;



    %% Do it
    rlist={};
    % Compute distance to the center of every atom
    for kk=1:size(digdata.Zdig,3)
        a_px = 2.68;
        R = digdata.Ratom{kk}-[digdata.Xc_px(kk);digdata.Yc_px(kk)];
        r = sqrt(R(1,:).^2+R(2,:).^2)/a_px;
        r=r(:);
        rlist{kk} = r;
    end
    minB =5;

    rMax = 110;                     % Max r   
    stepMax=floor(rMax/opts.BinStep);      % # of r steps
    redges = [0 minB];              % Small bin edge
    redges = [redges minB+[1:1:stepMax]*opts.BinStep]; % More bin edges
    rcen = (redges(1:end-1) + redges(2:end))/2;  

    nr = zeros(length(rcen),size(digdata.Zdig,3));

    areas = pi*redges.^2;areas=areas(:);
    dN = diff(areas);
    dN = dN(:);   

    for kk=1:size(digdata.Zdig,3)
        n = histcounts(rlist{kk},redges);
        nr(:,kk) = n(:)./dN;
    end

    nr_mean = mean(nr,2);
    nr_std = std(nr,0,2)/sqrt(size(nr,2));
    % pdf fit
    r_all = rlist{:};

    pdf_r = @(r,sigma) 2*pi*r.*exp(-r.^2/(2*sigma^2))/(2*pi*sigma^2);
    cdf_r = @(r,sigma) 1-exp(-r.^2/(2*sigma^2));
    N0=trapz([0 rcen],[mean(nr_mean(1:2)); nr_mean]);
    foo= @(s,r) sqrt(2/pi)*N0/s*exp(-r.^2/(2*s.^2));
    iMinGauss = find(flip(nr_mean)>0.1,1);

    if isempty(iMinGauss)
        rMinGauss=0;
    else
        rMinGauss = rcen(length(rcen)-iMinGauss);
    end
    rMinGauss=0;
    r_all_sub=r_all;
    r_all_sub(r_all_sub<=rMinGauss)=[];


    [pdf_r_vals,pdf_r_cints] = mle(r_all_sub,'pdf',pdf_r,'cdf',cdf_r,...
        'Start',[20], ...
    'LowerBound',0,'TruncationBounds',[rMinGauss inf]);
    sigmaR = pdf_r_vals(1);
    sigmaR_err =(pdf_r_cints(2,1)-pdf_r_cints(1,1))*0.5;


    sigmaR_um = sigmaR*0.527;
    sigmaR_err_um = sigmaR_err*0.527;
    T_HO = sigma2T_nK(sigmaR_um*1e-6);

    % Empircal CDF
    [f,x]=ecdf(r_all*0.527);

    %% String Stuff
    Z   = [digdata.Zdig];
    P   = [digdata.Params];

   
    strC = ['$(x_c,y_c):' ...
        '(' num2str(round(mean(digdata.Xc_px),0)) ',' ...
        num2str(round(mean(digdata.Yc_px),0)) ')~\mathrm{px},~'  ...
        '(' num2str(round(mean(digdata.Xc_site),0)) ',' ...
        num2str(round(mean(digdata.Yc_site),0)) ')~\mathrm{site}$'];

        strS = ['$(x_\sigma,y_\sigma):' ...
        '(' num2str(round(mean(digdata.Xs_um),1)) ',' ...
        num2str(round(mean(digdata.Ys_um),1)) ')~\mu\mathrm{m},~'  ...
        '(' num2str(round(mean(digdata.Xs_site),1)) ',' ...
        num2str(round(mean(digdata.Ys_site),1)) ')~\mathrm{site}$'];
    trap_str = ['$V(r)=0.5 m\omega^2r^2$' newline ...
        '$\omega=2\pi~' num2str(opts.TrapOmega/(2*pi)) '~\mathrm{Hz},t='...
        num2str(opts.Tunneling) '~\mathrm{Hz}$' newline ...
        '$\rightarrow r^*= ' ...
        num2str(round(rstar,1)) '~\mu\mathrm{m},T_\mathrm{HO}=' num2str(round(T_HO,1)) '~\mathrm{nK}$'];    
    %% Find Center   
    Xc  = [digdata.Xc_site];
    Yc  = [digdata.Yc_site];  
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
    % 
    %% Compute Average Image
    
    Npics = size(Z,3);    
    Zsub = zeros(length(ii(3):ii(4)),length(ii(1):ii(2)));

    for nn = 1:Npics
        Zsub(:,:,nn) = Z(ii(3):ii(4),ii(1):ii(2),nn);
    end    
    ZsubBar = mean(Zsub,3);
    % 
    % % Radial profile of average image
    % [rVec,radial_charge_mean,radial_charge_mean_dev,n]=radial_profile(ZsubBar,opts.BinStep);

    %% Compute Radial Profile of each Image
    % radial_charge = zeros(length(rVec),Npics);      %  Radial Charge Matrix
    % radial_charge_err = zeros(length(rVec),Npics);  %  Radial Charge  Deviation Matrix
    % 
    % % Iterate over each each
    % for ii = 1:Npics
    %     % Compute radial profile
    %     [rVec,charge,charge_dev,n]=...
    %         radial_profile(Zsub(:,:,ii),opts.BinStep);
    %     radial_charge(:,ii) = charge;
    %     radial_charge_err(:,ii) = charge_dev./sqrt(n);
    % end  
    % 
    % % Average Radial charge profile
    % radial_charge_mean      = mean(radial_charge,2);     
    % radial_charge_mean_std  = std(radial_charge,1,2);     
    % 

    %% Create radial potential vector
    
    % % Constants
    % h = 6.62607015*10^-34;
    % m = 39.96399848*1.66053906660*10^-27;
    % lam = (1054*10^-9*1054*10^-9*1064*10^-9)^(1/3);
    % aL = lam/2;
    % 
    % % 120mW XDT + 2.5ER request lattice trap frequencies - calibrated
    % % 02/29/24 (reanalyzed 04/11/24)
    % omega_x = 2*pi*67.3;
    % omega_y = 2*pi*60.2;
    % omega_bar = (omega_x*omega_y).^(1/2);
    % 
    % % Harmonic potential in Hz
    % PotentialVector = 0.5*m*omega_bar^2*aL^2.*rVec.^2/h;
    
    %% Assign to output
    
    % out = struct;
    % out.CroppedImages                  = Zsub;
    % out.AverageImage                   = ZsubBar;
    % out.SiteVector1                    = r(1):r(2);
    % out.SiteVector2                    = r(3):r(4);
    % out.RadialCenter                   = [Xcbar Ycbar];
    % out.RadialVector                   = rVec;
    % out.PotentialVector                = PotentialVector;
    % out.AverageOccupation              = radial_charge_mean;
    % out.AverageOccupationUncertainty   = radial_charge_mean_std/sqrt(Npics);

%% Plotting

    if ~isfield(opts,'Parent')
        opts.Parent = figure('color','w','Position',[300 100 1000 450],...
            'Name','Radial','NumberTitle','off');
        hF = opts.Parent;
    else
        hF = opts.Parent;
        for kk=1:length(hF.Children)
            delete(hF.Children(1))
        end
    end
    W = hF.Position(3);
    H = hF.Position(4);

    if isfield(opts,'FigLabel') && ~isempty(opts.FigLabel)
        tFig=uicontrol('style','text','string',opts.FigLabel,...
            'units','pixels','backgroundcolor',...
            'w','horizontalalignment','left','parent',hF);
        tFig.Position(4)=tFig.Extent(4);
        tFig.Position(3)=tFig.Extent(3);
        tFig.Position(1:2)=[1 1];
    end    

    
    % 2D Image of Average Image
    hF.UserData.Axes{1}=subplot(1,2,1,'parent',hF);
    ax1= hF.UserData.Axes{1};
    ax1.UserData.subplot_inds = [1 2 1];
    ax1.Units='pixels';
    [x0, y0, w, h]=getAxesPos(1,2,1,W,H);
    ax1.Position = [x0 y0 w h];
    
    if size(digdata.Zdig,3)>1
        imagesc(r(1):r(2),r(3):r(4),imboxfilt(ZsubBar,opts.BinStep))
        strBoxCar = ['boxcar avg: ' num2str(opts.BinStep) '\times' num2str(opts.BinStep)];
        text(.99,.01,strBoxCar,'units','normalized','fontsize',10,...
        'verticalalignment','top','horizontalalignment','right',...
        'color','r','parent',ax1);
        ax1.XAxisLocation='Top';
        set(ax1,'FontSize',8);
    else
        imagesc(r(1):r(2),r(3):r(4),ZsubBar)

    end
    xlabel(ax1,'position (sites)');
    ylabel(ax1,'position (sites)');    
    
    strSummary = ['$n=' num2str(Npics) ' $ images' newline ...
        '$N=' num2str(round(mean(digdata.Natoms))) '\pm' ...        
        num2str(round(std(digdata.Natoms,1))) '$ atoms' newline ...
        '$\sigma_r = ' num2str(round(sigmaR_um,1))  '\pm' ...
         num2str(round(sigmaR_err_um,1)) '\mu\mathrm{m}$'];

    text(1,1,strSummary,'units','pixels','fontsize',10,...
        'horizontalalignment','left','verticalalignment','bottom',...
        'color','r','interpreter','latex','parent',ax1);  
    colormap(ax1,'bone');
    axis(ax1,'equal');
    axis(ax1,'tight');

     text(.99,.99,[strC newline strS],'units','normalized','fontsize',12,...
        'verticalalignment','top','horizontalalignment','right',...
        'color','r','parent',ax1,'interpreter','latex','fontname','arial');
        ax1.XAxisLocation='Top';
        set(ax1,'FontSize',8);

    cc1=colorbar(ax1);
    cc1.Label.String = 'charge occupation';

    set(ax1,'ydir','normal','box','on','linewidth',1);
    
    if isfield(opts,'nMaxShow')
       caxis(ax1,[0 opts.nMaxShow]); 
    end
    % Average charge density
    hF.UserData.Axes{2}=subplot(2,2,2,'parent',hF);
    ax2=hF.UserData.Axes{2};
    ax2.Units='pixels';
    ax2.UserData.subplot_inds = [2 2 2];
    subplot_inds=ax2.UserData.subplot_inds;

    [x0, y0, w, h]=getAxesPos(subplot_inds(1),subplot_inds(2),subplot_inds(3),W,H);
    ax2.Position = [x0 y0 w h];

    errorbar(rcen,nr_mean,nr_std,'ko','linewidth',1,'markerfacecolor',[.5 .5 .5],...
        'parent',ax2);
     hold(ax2,'on')
    xlabel(ax2,'radial distance (sites)');
    ylabel(ax2,'$n(r)$','interpreter','latex','fontsize',14);
    text(.99,.99,strSummary,'units','normalized','fontsize',12,'horizontalalignment','right',...
        'verticalalignment','top','interpreter','latex','parent',ax2);
    
    if isfield(opts,'nMaxShow')
       ylim(ax2,[0 opts.nMaxShow]); 
    end
    
    if isfield(opts,'rMaxShow')
       xlim(ax2,[0 opts.rMaxShow]); 
    end
    
    strRadialBin = ['radial bin ' char(916) 'r:' num2str(opts.BinStep)];
    text(.99,.01,strRadialBin,'units','normalized','horizontalalignment','right',...
        'verticalalignment','bottom','fontsize',10,'parent',ax2);
    rVec=linspace(0,100,100);
    plot(rVec,foo(sigmaR,rVec),'-')
    set(ax2,'box','on','linewidth',1,'fontsize',12,...
        'yaxislocation','right')

    % Plot radial average ndet std
    hF.UserData.Axes{3}=subplot(2,2,4,'parent',hF);
    ax3=hF.UserData.Axes{3};
    ax3.Units='pixels';
    ax3.UserData.subplot_inds = [2 2 4];
    [x0, y0, w, h]=getAxesPos(2,2,4,W,H);
    ax3.Position = [x0 y0 w h];
    c=[0    0.4470    0.7410];

%     errorbar(1e-3*site2Vr_Hz(rcen),nr_mean,nr_std,'o-','linewidth',1, ...
%         'markerfacecolor',c,'markeredgecolor',c*.5,'linewidth',1,'parent',ax3)
%     xlim(ax3,[0 10]);
%     xlabel(ax3,'$V(r)~(\mathrm{kHz})$','interpreter','latex','fontsize',14)
%     ylabel(ax3,'$n(r)$','interpreter','latex','fontsize',14);
%     set(ax3,'box','on','linewidth',1,'fontsize',12,...
%         'yaxislocation','right')

    plot(x,100*f,'k-','linewidth',2,'parent',ax3);
    xlabel(ax3,'radial position (\mum)')
    xlim(ax3,[0 max(x)]);
    ylim([0 102]);
    set(ax3,'box','on','linewidth',1,'fontsize',12,...
        'yaxislocation','right')
    ylabel('% within radius')
    hold(ax3,'on');
        

    plot([1 1]*rstar,[0 1]*interp1(x(2:end),100*f(2:end),rstar),'k--','linewidth',1);
    plot([rstar x(end)],[1 1]*interp1(x(2:end),100*f(2:end),rstar),'k--','linewidth',1);
text(.99,.01,trap_str,'interpreter','latex','horizontalalignment','right',...
    'verticalalignment','bottom','fontsize',12,'units','normalized')
    % 
    % plot(rVec(2:end),radial_charge_mean_std(2:end)/sqrt(Npics),...
    %     'ko','markerfacecolor',[.5 .5 .5],...
    %     'markersize',8,'linewidth',1,'parent',ax3);
    % title(ax3,'standard error');
    % ylabel(ax3,'std $\{N(r)\}/\sqrt{n}$','interpreter','latex','fontsize',14);
    % xlabel(ax3,'radial distance (sites)');    
    % text(.99,.99,strSummary,'units','normalized','fontsize',12,'horizontalalignment','right',...
    %     'verticalalignment','top','interpreter','latex','parent',ax3);
    % 
    % if isfield(opts,'rMaxShow')
    %    xlim(ax3,[0 opts.rMaxShow]); 
    % end
    % text(.01,.99,strRadialBin,'units','normalized','horizontalalignment','left',...
    %     'verticalalignment','top','fontsize',10,'parent',ax3);
    % ax3.YAxisLocation='right';


    hF.UserData.Axes{1}.UserData.subplot_inds = [1 2 1];
    hF.UserData.Axes{2}.UserData.subplot_inds = [2 2 2];
    hF.UserData.Axes{3}.UserData.subplot_inds = [2 2 4];

    hF.SizeChangedFcn=@figResize;
  
end
  function figResize(src,evt)
        for kk=1:length(src.UserData.Axes)
            ax=src.UserData.Axes{kk};
            inds=ax.UserData.subplot_inds;
            [a b c d]=getAxesPos(inds(1),inds(2),inds(3),src.Position(3),src.Position(4));
            src.UserData.Axes{kk}.Position=[a b c d];                
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


function [axX,axY,axWidth,axHeight]=getAxesPos(num_rows,num_cols,index,figW,figH)
yTop=60;
yBot=50;

xLeft=60;
xRight=65;

ySpace=60;
xSpace=60;

rowNumber = floor((index-1)/num_cols)+1;
colNumber = mod(index-1,num_cols)+1;

axHeight = (figH-yTop-yBot-ySpace*(num_rows-1))/num_rows;
axWidth = (figW-xLeft-xRight-xSpace*(num_cols-1))/num_cols;

axX = xLeft + (axWidth+xSpace)*(colNumber-1);
axY=(figH-yTop-axHeight)-(rowNumber-1)*(axHeight+ySpace);

end

