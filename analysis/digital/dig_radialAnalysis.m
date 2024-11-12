function [hF,output] = dig_radialAnalysis(digdata,opts)
if nargin == 1;opts = struct;end    

% How to calculate n(r)
if ~isfield(opts,'BinStep');        opts.BinStep = 3;end   
if ~isfield(opts,'Bin0');           opts.Bin0 = 5;end   
if ~isfield(opts,'rMax');           opts.rMax = 110;end

% Fittng n(r)
if ~isfield(opts,'GaussFitDensityMax');opts.GaussFitDensityMax = [.05 1];end
%% Parameters Relvant do Fermi-Hubbard Model

if ~isfield(opts,'TrapOmega');      opts.TrapOmega = [];end
if ~isfield(opts,'Tunneling');      opts.Tunneling = [];end
if ~isfield(opts,'U');              opts.Tunneling = [];end


if ~isfield(opts,'TrapOmega');      opts.TrapOmega   = 2*pi*65;end
if ~isfield(opts,'Tunneling');      opts.Tunneling   = 563;end
if ~isfield(opts,'Interaction');    opts.Interaction = 1e3;end

%% Constants

% Physical constants
h = 6.62607015e-34;         % planck's constant [Js]
amu = 1.66053906660e-27;    % atomic mass unit[kg]    
m = amu*39.96399848;        % 40K mass [kg]
aL = (1054e-9)/2;           % lattice spacing [m]
aL_um = aL*1e6;
kB= 1.380649e-23;           % J/K

%% Trap Potential V(r) functions

site2Vr_Hz =        @(n) 0.5*m*opts.TrapOmega^2*aL^2*n.^2/h;
sigma2T_nK =        @(sigmaR) 1e9*sigmaR^2*m*opts.TrapOmega^2/kB;
sigma2Tovert =      @(sigmaR) sigmaR^2*m*opts.TrapOmega^2/(h*opts.Tunneling);

%% Tunnel Connection Radius

nstar = linspace(0,100,1e4);
istar=find([site2Vr_Hz(nstar)]>=4*opts.Tunneling,1);
rstar=nstar(istar)*aL_um;

if isfield(digdata,'Lattice_px')
    a_px = digdata.Lattice_px;
else
    a_px = 2.68;
end

%% Get Data
Z   = [digdata.Zdig];
Zbar = mean(Z,3);
P   = [digdata.Params];

%% Calculate radial profile
rlist={};
% Compute distance to the center of every atom
for kk=1:size(digdata.Zdig,3)
    R = digdata.Ratom{kk}-[digdata.Xc_px(kk);digdata.Yc_px(kk)];
    r = sqrt(R(1,:).^2+R(2,:).^2)/a_px;
    r=r(:);
    rlist{kk} = r;
end

% Calculate the bins
% Around r=0, find all a points within pi*Bin0^2
% For all other radii bin are separated by BinStep

stepMax=floor(opts.rMax/opts.BinStep);                   % # of r steps
redges = [0 opts.Bin0];                                  % About r=0
redges = [redges opts.Bin0+[1:1:stepMax]*opts.BinStep];  % Rest of edges
rcen = (redges(1:end-1) + redges(2:end))/2;             % bin centers

% Jacobian via ideal circle
areas = pi*redges.^2;areas=areas(:);                % Area of each bin edge
dN = diff(areas);                                   % Jacobian
dN = dN(:);                                         % Resize

% Jacobian via numerical circle
n = [-opts.rMax:1:opts.rMax];
[nn1,nn2]=meshgrid(n,n);
R=(nn1.^2+nn2.^2).^(1/2);
Nr = zeros(length(redges),1);
for rr=1:length(redges)
    Nr(rr)=sum(R(:)<redges(rr));
end
dN2 = diff(Nr);

% Calculate the histograms
nr = zeros(length(rcen),size(digdata.Zdig,3));
for kk=1:size(digdata.Zdig,3)
    n = histcounts(rlist{kk},redges);
    nr(:,kk) = n(:)./dN2;       % normalize by discrete jacboian
    % nr(:,kk) = n(:)./dN;      % normalize by continuous jacboian
end

% take the average
nr_mean = mean(nr,2);
nr_std = std(nr,0,2)/sqrt(size(nr,2));

%% Get all r for all atoms
r_all = rlist{:};               % r for all atoms
[f,x]=ecdf(r_all*aL_um);        % Empircal CDF
%% Calculate numerical second moment
% The gaussian radius assuming a symmetric gaussian distribution is
% related to the expectation value of the radius.

% integral (2*pi*r)*P(r)*r = sigma_R*sqrt(pi/2)
r_expect = mean(r_all,'all');            % expected r
sigma_r_numerical = r_expect*sqrt(2/pi); % expected r to sigma_R

T_HOt = sigma2Tovert(sigma_r_numerical*aL);
T_HO_nK = 1e9*sigma2Tovert(sigma_r_numerical*aL)*h*opts.Tunneling/kB;

%% Fit Data to Gaussian PDf
% Fit the measured atoms to a probability density function
rMin = zeros(length(opts.GaussFitDensityMax),1);
sigma_r_gauss_fit = zeros(length(opts.GaussFitDensityMax),1);
sigma_r_gauss_fit_err = zeros(length(opts.GaussFitDensityMax),1);

% Numerical integral of n(r)
I0=trapz([0 rcen],[mean(nr_mean(1:2)); nr_mean]);

% Iterate over max density fit
for nn=1:length(opts.GaussFitDensityMax)
    iMinFlip = find(flip(nr_mean)>opts.GaussFitDensityMax(nn),1);    
    if isempty(iMinFlip)
        iMin = 1;           % 
        rMin(nn)=0;         % Fit all radii
        Iouter(nn) = I0;    % Total integral is the max
    else
        iMin = length(rcen)-iMinFlip; % index at which n(r(ind))=n_limit
        rMin(nn) = rcen(iMin);% radies at which n(rMin)=n_limit
        Iouter(nn) = trapz(rcen(iMin:end),nr_mean(iMin:end)); % partial integral
    end
    % Fit the distribution
    [sigma_r_gauss_fit(nn),sigma_r_gauss_fit_err(nn)]=pdf_gauss_fit(r,rMin(nn));
end

% Convert the fit into n(r)
nr_partial_cdf  = @(s,R) 1-erf(R./sqrt(2*s^2)); % integral of foo [R,infinity]
nr_fit          = @(s,r,N0) sqrt(2/pi)/s*exp(-r.^2/(2*s.^2))*N0;
     

%% String Stuff

% Center of the cloud
strC = ['$(x_c,y_c):' ...
    '(' num2str(round(mean(digdata.Xc_px),0)) ',' ...
    num2str(round(mean(digdata.Yc_px),0)) ')~\mathrm{px},~'  ...
    '(' num2str(round(mean(digdata.Xc_site),0)) ',' ...
    num2str(round(mean(digdata.Yc_site),0)) ')~\mathrm{site}$'];

% Second moment of cloud
strS = ['$(x_\sigma,y_\sigma):' ...
    '(' num2str(round(mean(digdata.Xs_um),1)) ',' ...
    num2str(round(mean(digdata.Ys_um),1)) ')~\mu\mathrm{m},~'  ...
    '(' num2str(round(mean(digdata.Xs_site),1)) ',' ...
    num2str(round(mean(digdata.Ys_site),1)) ')~\mathrm{site}$'];
    

    if ~isnan(opts.TrapOmega) && ~isnan(opts.Tunneling)
        trap_str = ['$\omega=2\pi\cdot' num2str(opts.TrapOmega/(2*pi)) '~\mathrm{Hz},t='...
            num2str(opts.Tunneling) '~\mathrm{Hz}$'];    
        trap_str = [trap_str newline ...
            '$\rightarrow r^*= ' ...
            num2str(round(rstar,1)) '~\mu\mathrm{m},' ...
            'T\approx m\omega^2\langle\sigma_r\rangle^2=' num2str(T_HOt,'%.1f') ...
           't~(' num2str(round(T_HO_nK)) '~\mathrm{nK})$'];    
    else
        trap_str=[];
    end

    lineN = ['$' num2str(round(mean(digdata.Natoms))) '\pm' ...        
        num2str(round(std(digdata.Natoms,1))) '$ atoms (n=' num2str(size(Z,3)) ')'];
    lineSigmaExpect = ['$\langle \sigma_r \rangle=' ...
        num2str(round(aL_um*sigma_r_numerical,1)) '~\mu\mathrm{m}$'];

    strRadialSummary = [ lineN newline ...
        lineSigmaExpect];

    %% Find Center and limits in site space
    n1 = digdata.n1;
    n2 = digdata.n2;        
    [nn1,nn2]=meshgrid(n1,n2);

    for nn=1:size(Z,3)
        Zthis = Z(:,:,nn);
        n1c(nn) = sum(Zthis.*nn1,'all')/sum(Zthis,'all');
        n2c(nn) = sum(Zthis.*nn2,'all')/sum(Zthis,'all');    
        n1c_sq(nn) = sum(Zthis.*nn1.^2,'all')/sum(Zthis,'all');
        n2c_sq(nn) = sum(Zthis.*nn2.^2,'all')/sum(Zthis,'all');    
        n1_sigma(nn) = sqrt(n1c_sq(nn)-n1c(nn)^2);
        n2_sigma(nn) = sqrt(n2c_sq(nn)-n2c(nn)^2);
    end
    n_sigma_med =[median(n1_sigma) median(n2_sigma)];
    nc_med = [median(n1c) median(n2c)];    
    n1_lim = nc_med(1)+4.5*[-1 1]*n_sigma_med(1);
    n2_lim = nc_med(2)+4.5*[-1 1]*n_sigma_med(2);

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


    %% Assign to output
    output = struct;
    output.Z                        = Zbar;
    output.NumImages                = size(Z,3);
    output.RadialOccupation         = nr;
    output.RadialOccupationMean     = nr_mean;
    output.RadialOccupationError    = nr_std;

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
        

    if isfield(digdata,'SourceDirectory') && ...
            ~isempty(digdata.SourceDirectory) &&  ...
            ~isequal(digdata.SourceDirectory,{'GUI'})
        tstr = digdata.SourceDirectory;
    else
        tstr  = digdata.FileNames{1};
    end    
    
        tFig=uicontrol('style','text','string',tstr,...
            'units','pixels','backgroundcolor',...
            'w','horizontalalignment','left','parent',hF);
        tFig.Position(4)=tFig.Extent(4);
        tFig.Position(3)=tFig.Extent(3);
        tFig.Position(1:2)=[1 1];
    
    % 2D Image of Average Image
    hF.UserData.Axes{1}=subplot(1,2,1,'parent',hF);
    ax1= hF.UserData.Axes{1};
    ax1.UserData.subplot_inds = [1 2 1];
    ax1.Units='pixels';
    [x0, y0, w, h]=getAxesPos(1,2,1,W,H);
    ax1.Position = [x0 y0 w h];
    
    if size(digdata.Zdig,3)>1
        imagesc(n1,n2,imboxfilt(Zbar))
        strBoxCar = ['boxcar avg: ' num2str(opts.BinStep) '\times' num2str(opts.BinStep)];
        text(.99,.01,strBoxCar,'units','normalized','fontsize',10,...
        'verticalalignment','top','horizontalalignment','right',...
        'color','r','parent',ax1);
        ax1.XAxisLocation='Top';
        set(ax1,'FontSize',8);
    else
        imagesc(n1,n2,Zbar)
    end
    xlabel(ax1,'position (sites)');
    ylabel(ax1,'position (sites)');    



    text(1,1,strRadialSummary,'units','pixels','fontsize',10,...
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
    xlim(ax1,n1_lim);
    ylim(ax1,n2_lim);

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
    text(.99,.99,strRadialSummary,'units','normalized','fontsize',12,'horizontalalignment','right',...
        'verticalalignment','top','interpreter','latex','parent',ax2);
    
    if isfield(opts,'nMaxShow')
       ylim(ax2,[0 opts.nMaxShow]); 
    end
    
    if isfield(opts,'rMaxShow')
       xlim(ax2,[0 opts.rMaxShow]); 
    end
    
    strRadialBin = ['radial bin ' char(916) 'r:' num2str(opts.BinStep)];
    text(.01,1,strRadialBin,'units','normalized','horizontalalignment','left',...
        'verticalalignment','bottom','fontsize',8,'parent',ax2);

    rVec=linspace(0,100,100);
    for nn=1:length(opts.GaussFitDensityMax)
        s=sigma_r_gauss_fit(nn);
        s_err=sigma_r_gauss_fit_err(nn);

        N = Iouter(nn)/nr_partial_cdf(s,rMin(nn));
        pGaussFits(nn)=plot(rVec,nr_fit(s,rVec,N),'-');
        legStr{nn}=['$' num2str(s*aL_um,'%.1f')  ...
            '(' num2str(round(10*s_err*aL_um)) ')~\mu\mathrm{m}$' ...
            ' $n<' num2str(opts.GaussFitDensityMax(nn)) '$'];
    end
    if ~isempty(opts.GaussFitDensityMax)
    legend(pGaussFits,legStr,'interpreter','latex','fontsize',8,...
        'location','southeast');
    end

    set(ax2,'box','on','linewidth',1,'fontsize',12,...
        'yaxislocation','right')

    % Plot radial average ndet std
    hF.UserData.Axes{3}=subplot(2,2,4,'parent',hF);
    ax3=hF.UserData.Axes{3};
    ax3.Units='pixels';
    ax3.UserData.subplot_inds = [2 2 4];
    [x0, y0, w, h]=getAxesPos(2,2,4,W,H);
    ax3.Position = [x0 y0 w h];
    
    plot(x,100*f,'k-','linewidth',2,'parent',ax3);
    xlabel(ax3,'radial distance (\mum)')
    xlim(ax3,[0 max(x)]);
    ylim([0 102]);
    set(ax3,'box','on','linewidth',1,'fontsize',12,...
        'yaxislocation','right')
    ylabel('% within radius')
    hold(ax3,'on');
        

    plot([1 1]*rstar,[0 1]*interp1(x(2:end),100*f(2:end),rstar),'k--','linewidth',1);
    plot([rstar x(end)],[1 1]*interp1(x(2:end),100*f(2:end),rstar),'k--','linewidth',1);

    if ~isempty(trap_str)    
        text(.98,.02,trap_str,'interpreter','latex','horizontalalignment','right',...
        'verticalalignment','bottom','fontsize',10,'units','normalized',...
        'backgroundcolor','w','margin',1)
    end
    hF.UserData.Axes{1}.UserData.subplot_inds = [1 2 1];
    hF.UserData.Axes{2}.UserData.subplot_inds = [2 2 2];
    hF.UserData.Axes{3}.UserData.subplot_inds = [2 2 4];

    hF.SizeChangedFcn=@figResize;
  
end

function [sigmaR,sigmaR_err]=pdf_gauss_fit(r,rMin)

    % Probability density function
    pdf_r = @(r,sigma) 2*pi*r.*exp(-r.^2./(2*sigma.^2))./(2*pi*sigma.^2);
    % Cummulate density function (goes from 0 to 1);
    cdf_r = @(r,sigma) 1-exp(-r.^2/(2*sigma^2));


    r_expect = mean(r,'all');   % expected r
    r(r<rMin)=[];                   % remove data points smaller rMin

    % Fit it
    [pdf_r_vals,pdf_r_cints] = mle(r,...
        'pdf',pdf_r,'cdf',cdf_r,...
        'Start',r_expect*sqrt(2/pi), ...
        'LowerBound',0,'TruncationBounds',[rMin inf]);
    sigmaR = pdf_r_vals(1);
    sigmaR_err =(pdf_r_cints(2,1)-pdf_r_cints(1,1))*0.5;
end


  function figResize(src,evt)
        for kk=1:length(src.UserData.Axes)
            ax=src.UserData.Axes{kk};
            inds=ax.UserData.subplot_inds;
            [a b c d]=getAxesPos(inds(1),inds(2),inds(3),src.Position(3),src.Position(4));
            src.UserData.Axes{kk}.Position=[a b c d];                
        end
  end
  
  % Obsolete
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

