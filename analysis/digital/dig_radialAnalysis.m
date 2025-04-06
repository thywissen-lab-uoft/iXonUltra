function [hF,output] = dig_radialAnalysis(digdata,opts)
% Author : CJ Fujiwara
%
% This code takes digitized data and does various radial analysis which
% assumes cylindrical symmetry, which is quite common.
if nargin == 1;opts = struct;end    

hubbard = load('hubbard.mat');
hubbard=hubbard.hubbard;

v = hubbard.depth;
t11 = hubbard.fr*reshape(hubbard.Tunneling(1,1,:),1,[]);
t12 = hubbard.fr*reshape(hubbard.Tunneling(1,2,:),1,[]);
t13 = hubbard.fr*reshape(hubbard.Tunneling(1,3,:),1,[]);
t21 = hubbard.fr*reshape(hubbard.Tunneling(2,1,:),1,[]);
t22 = hubbard.fr*reshape(hubbard.Tunneling(2,2,:),1,[]);
t23 = hubbard.fr*reshape(hubbard.Tunneling(2,3,:),1,[]);

%% Default Binning Options
% These are the default bin settings.

if ~isfield(opts,'BinStep');        opts.BinStep = 4;end   
if ~isfield(opts,'Bin0');           opts.Bin0 = 4;end   
if ~isfield(opts,'rMax');           opts.rMax = 110;end

strRadialBin = ['radial bin ' char(916) ...
    'r:' num2str(opts.BinStep) ', R_0: ' num2str(opts.Bin0)];

% YOU NEED A BETTER VERSION OF MATLAB TO DO OTHER DENSITY MAXIMUMS
if ~isfield(opts,'GaussFitDensityMax');opts.GaussFitDensityMax = [1];end
%% Constants

% Physical constants
h       = 6.62607015e-34;         % planck's constant [Js]
amu     = 1.66053906660e-27;      % atomic mass unit[kg]    
m       = amu*39.96399848;        % 40K mass [kg]
aL      = (1054e-9)/2;            % lattice spacing [m]
aL_um   = aL*1e6;                 % lattice spacign [m]
kB      = 1.380649e-23;           % boltzmann constant [J/K]

%% Parameters Relvant do Fermi-Hubbard Model
% Default Tunneling Values


if ~isfield(opts,'TrapOmega');      opts.TrapOmega   = 2*pi*67;end%03/05 2*pi*52*sqrt(1.0965);end % from 11/25/24
if ~isfield(opts,'Tunneling');      opts.Tunneling   = 563;end
if ~isfield(opts,'Interaction');    opts.Interaction = 1e3;end


%opts.TrapOmega = 2*pi*90;
%opts.Tunnelling = 435;

input_data{1,1} = ['V₀ [Eᵣ]'];
input_data{2,1} = ['ωr/(2π) [Hz]'];
input_data{3,1} = ['ωz/(2π) [Hz]'];
input_data{4,1} = ['B [Gauss]'];

input_data{1,2} = 2.5;
input_data{2,2} = 67;
input_data{3,2} = 67*5;
input_data{4,2} = 195;

output_data{1,1} = 'tunneling [Hz]';
output_data{1,2} = '[563, 20, .1]';
output_data{2,1} = 'band gap [Hz]';
output_data{2,2} = '[2000, 200, 100]';
output_data{3,1} = 'tunnel r [sites]';
output_data{3,2} = '20';
output_data{4,1} = 'tunnel z [sites]';
output_data{4,2} = '5';
output_data{5,1} = 'U/t';
output_data{5,2} = '5';

    function output = calcHubbard
        % v0 = params.depth;
        v0=2.5;
        omega_r = 2*pi*67;
        omega_z = 2*pi*67*5;

        % First Band Tunneling
        tunnel_11 = interp1(v,t11,v0);  % one-site
        tunnel_12 = interp1(v,t12,v0);  % two-site
        tunnel_13 = interp1(v,t13,v0);  % three-site
        output{1,1}='s-tunneling [Hz]';
        output{1,2}=['(' num2str(round(tunnel_11,0)) ',' ...
            num2str(round(tunnel_12,0)) ',' ...
            num2str(round(tunnel_13,0)) ')'];

        % Second Band Tunneling
        % tunnel_21 = interp1(v,t21,v0);  % one-site
        % tunnel_22 = interp1(v,t22,v0);  % two-site
        % tunnel_23 = interp1(v,t23,v0);  % three-site
        % output{2,1}='p-tunneling [Hz]';
        % output{2,2}=['(' num2str(round(tunnel_21,0)) ',' ...
        %     num2str(round(tunnel_22,0)) ',' ...
        %     num2str(round(tunnel_23,0)) ')'];

        % Band Gap
        bg_1d = interp1(v,hubbard.fr*hubbard.BandGap1D,v0);
        bg_2d = interp1(v,hubbard.fr*hubbard.BandGap2D,v0);
        bg_3d = interp1(v,hubbard.fr*hubbard.BandGap3D,v0);
        output{2,1}='1D,2D,3D Gap [Hz]';
        output{2,2}=['(' num2str(round(bg_1d,0)) ',' ...
            num2str(round(bg_2d,0)) ',' ...
            num2str(round(bg_3d,0)) ')'];

        % Tunneling Radius
        r_tunnel = sqrt(4*tunnel_11*h/(0.5*m*omega_r^2));
        z_tunnel = sqrt(4*tunnel_11*h/(0.5*m*omega_z^2));
        output{3,1} = 'r tunnel [sites]';
        output{3,2} = num2str(round(r_tunnel/aL,1));
        output{4,1} = 'z tunnel [sites]';
        output{4,2} = num2str(round(z_tunnel/aL,1));

        % Hubbard U
        output{5,1} = 'U/t';
        output{5,2} = 'TBD';
    end
output_data = calcHubbard();

%% Trap Potential V(r) functions

site2Vr_Hz =        @(n) 0.5*m*opts.TrapOmega^2*aL^2*n.^2/h;
sigma2T_nK =        @(sigmaR) 1e9*sigmaR^2*m*opts.TrapOmega^2/kB;
sigma2Tovert =      @(sigmaR) sigmaR^2*m*opts.TrapOmega^2/(h*opts.Tunneling);

%% Tunnel Connection Radius

nstar   = linspace(0,100,1e4);
istar   = find([site2Vr_Hz(nstar)]>=4*opts.Tunneling,1);
rstar   = nstar(istar)*aL_um;

if isfield(digdata,'Lattice_px')
    a_px = digdata.Lattice_px;
else
    a_px = 2.68;
end

%% Calculate Radial Jacobian and Bins

stepMax     = floor(opts.rMax/opts.BinStep);                  % # of r steps
redges      = [0 opts.Bin0];                                  % About r=0
redges      = [redges opts.Bin0+[1:1:stepMax]*opts.BinStep];  % Rest of edges
rcen        = (redges(1:end-1) + redges(2:end))/2;            % bin centers

% Jacobian via ideal circle
areas           = pi*redges.^2;areas=areas(:);                % Area of each bin edge
jacobian_circle = diff(areas);                                % Jacobian
jacobian_circle = jacobian_circle(:);                         % Resize

% Jacobian via numerical circle
n               = [-opts.rMax:1:opts.rMax];
[nn1,nn2]       = meshgrid(n,n);
R               = (nn1.^2+nn2.^2).^(1/2);
Nr              = zeros(length(redges),1);
r_center        = zeros(length(redges),1);

for rr=1:length(redges)
    indsLess    = [R<redges(rr)];    
    Nr(rr)      = sum(indsLess,'all');
    if rr==1
        indsMore    = [R>0];
    else
        indsMore    = [R>=redges(rr-1)];
    end
    indsSlice    = logical(indsLess.*indsMore);
    r_center(rr) = mean(R(indsSlice));  
end
r_center        = r_center(2:end)';
rcen            = r_center;
jacobian_numerical = diff(Nr);


%% Get Data
Z           = [digdata.Zdig];
Z           = Z(:,:,:,1);
Zbar        = mean(Z,3);
P           = [digdata.Params];

%% Calculate radial profile
% For every image compute distance to the center of every atom
rlist={};
img_ind = 1; % Image index (in case of double shot, we always use first)
for kk=1:size(digdata.Zdig,3)
    Rc =  [digdata.Xc_px(kk,img_ind);digdata.Yc_px(kk,img_ind)];
    R = digdata.Ratom{kk,img_ind}-Rc;
    r = sqrt(R(1,:).^2+R(2,:).^2)/a_px;
    r=r(:);
    rlist{kk} = r;
end
r_all = rlist{:};               % r for all atoms
[f,x] = ecdf(r_all*aL_um);        % Empircal CDF (cummulative density function)
%% Calculate Radial Histogram
nr = zeros(length(rcen),size(digdata.Zdig,3));
for kk=1:size(digdata.Zdig,3)
    n = histcounts(rlist{kk},redges);
    nr(:,kk) = n(:)./jacobian_numerical;       % normalize by discrete jacboian
    % nr(:,kk) = n(:)./jacobian_circle;      % normalize by continuous jacboian
end

% Average over each image
nr_mean = mean(nr,2);
nr_std = std(nr,0,2)/sqrt(size(nr,2));

%% Fit Observations to 2D Gaussian PDF
% Fit the measured atoms to a 2D Gaussian probability density function
rMin                    = zeros(length(opts.GaussFitDensityMax),1);
gauss_sigma             = zeros(length(opts.GaussFitDensityMax),2);

% Numerical integral of n(r)
I0=trapz([0 rcen],[mean(nr_mean(1:2)); nr_mean]);

% Iterate over max density fit
gauss_str={};
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
    [gauss_sigma(nn,1),gauss_sigma(nn,2)]=pdf_gauss_fit(r_all,rMin(nn));
    s=gauss_sigma(nn,1);

    T_HO = m*(2*pi*67)^2*(gauss_sigma(nn,1)*aL)^2/kB;
    T_HO_t = (kB*T_HO)/(h*563);
    T_HO_nK = T_HO*1e9;
    gauss_str{nn}=['Gauss : $\sigma=' num2str(round(s),'%.1f')  ...
         '~\mathrm{sites}~(' num2str(round(s*aL_um,1)) '~\mu \mathrm{m})' ...
         ';T=' num2str(round(T_HO_nK,0)) '~\mathrm{nK}~(' num2str(round(T_HO_t,1)) 't)'...
        '~[\mathrm{fit}~n(r)<' num2str(opts.GaussFitDensityMax(nn)) ']$'];
end

% Convert the fit into n(r)
nr_partial_cdf  = @(s,R) 1-erf(R./sqrt(2*s^2)); % integral of foo [R,infinity]
nr_fit          = @(s,r,N0) sqrt(2/pi)/s*exp(-r.^2/(2*s.^2))*N0;

%% Calculate numerical second moment
% The gaussian radius assuming a symmetric gaussian distribution is
% related to the expectation value of the radius.

% integral (2*pi*r)*P(r)*r = sigma_R*sqrt(pi/2)
r_expect = mean(r_all,'all');            % expected r
sigma_r_numerical = r_expect*sqrt(2/pi); % expected r to sigma_R

T_HOt = sigma2Tovert(sigma_r_numerical*aL);
T_HO_nK = 1e9*sigma2Tovert(sigma_r_numerical*aL)*h*opts.Tunneling/kB;
%% Fit Data to Gibb's estimate
% CJF : I am not familiar with this fitting protocol as it is something
% that RL and JHT came up with. This should really be renamed to something
% that is more clear into the physical assumptions. Example "Gibbs" I don't
% think really means anything besides the fact it uses a Gibbs ensemble of
% some kind....
%
% Also a note.. is that this function is a FIT function, and not a
% probability density fit. So this will depend heavily on binning. It is
% more correct to fit the probability. Please fix this.

% Guesses
Gz0 = 1;
Gs = sigma_r_numerical;
w       = 1./(nr_std.^2); % weights as 1/variance
w(isinf(w)|isnan(w)) = 0;


gibbsFit=fittype('2./(z0*exp(-r.^2./(2*s.^2)) + z0.^(-1)*exp(r.^2./(2*s.^2)) + 2)','independent',{'r'},...
    'coefficients',{'z0','s'});
options=fitoptions(gibbsFit);          
options.TolFun      = 1E-14;
options.Lower       = [-10 0];
options.Upper       = [30 2*Gs];            
options.StartPoint  = [Gz0 Gs];     
options.MaxIter     = 3000;
options.MaxFunEvals = 3000;

        
GibbsFit=fit(rcen',nr_mean,gibbsFit,options);
cint = confint(GibbsFit,0.683);
Gz0Err = (cint(2,1)-cint(1,1))*0.5;
GsErr = (cint(2,2)-cint(1,2))*0.5;

T_HOt_g = sigma2Tovert(GibbsFit.s*aL);
T_HO_g_nK = 1e9*sigma2Tovert(GibbsFit.s*aL)*h*opts.Tunneling/kB;

npeak_g = feval(GibbsFit,0);
n_doublon = (GibbsFit.z0).^2./((GibbsFit.z0).^2 + 2*(GibbsFit.z0) + 1);
n_up = npeak_g/2 + n_doublon;
     

%% String Stuff

% Center of the cloud
strC = ['$(x_c,y_c):' ...
    '(' num2str(round(mean(digdata.Xc_px(:,1)),0)) ',' ...
    num2str(round(mean(digdata.Yc_px(:,1)),0)) ')~\mathrm{px},~'  ...
    '(' num2str(round(mean(digdata.Xc_site(:,1)),0)) ',' ...
    num2str(round(mean(digdata.Yc_site(:,1)),0)) ')~\mathrm{site}$'];

% Second moment of cloud
strS = ['$(x_\sigma,y_\sigma):' ...
    '(' num2str(round(mean(digdata.Xs_um(:,1)),1)) ',' ...
    num2str(round(mean(digdata.Ys_um(:,1)),1)) ')~\mu\mathrm{m},~'  ...
    '(' num2str(round(mean(digdata.Xs_site(:,1)),1)) ',' ...
    num2str(round(mean(digdata.Ys_site(:,1)),1)) ')~\mathrm{site}$'];


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

lineN = ['$' num2str(round(mean(digdata.Natoms(:,1)))) '\pm' ...        
    num2str(round(std(digdata.Natoms(:,1),1))) '$ atoms (n=' num2str(size(Z,3)) ')'];
lineSigmaExpect = ['$\langle \sigma_r \rangle=' ...
    num2str(round(aL_um*sigma_r_numerical,1)) '~\mu\mathrm{m}$'];

strRadialSummary = [ lineN newline ...
    lineSigmaExpect];

%% Find Center and limits in site space
n1 = digdata.n1;
n2 = digdata.n2;        
[nn1,nn2]=meshgrid(n1,n2);

for nn=1:size(Z,3)
    Zthis           = Z(:,:,nn);
    n1c(nn)         = sum(Zthis.*nn1,'all')/sum(Zthis,'all');
    n2c(nn)         = sum(Zthis.*nn2,'all')/sum(Zthis,'all');    
    n1c_sq(nn)      = sum(Zthis.*nn1.^2,'all')/sum(Zthis,'all');
    n2c_sq(nn)      = sum(Zthis.*nn2.^2,'all')/sum(Zthis,'all');    
    n1_sigma(nn)    = sqrt(n1c_sq(nn)-n1c(nn)^2);
    n2_sigma(nn)    = sqrt(n2c_sq(nn)-n2c(nn)^2);
end
n_sigma_med =[median(n1_sigma) median(n2_sigma)];
nc_med = [median(n1c) median(n2c)];    
n1_lim = nc_med(1)+3*[-1 1]*n_sigma_med(1)+[-5 5];
n2_lim = nc_med(2)+3*[-1 1]*n_sigma_med(2)+[-5 5];

n1_lim(1) = max([n1_lim(1) n1(1)]);
n1_lim(2) = min([n1_lim(2) n1(end)]);

n2_lim(1) = max([n2_lim(1) n2(1)]);
n2_lim(2) = min([n2_lim(2) n2(end)]);


%% Assign to output
output = struct;
output.Z                        = Zbar;
output.NumImages                = size(Z,3);
output.RadialVec                = rcen'.*aL_um;
output.RadialOccupation         = nr;
output.RadialOccupationMean     = nr_mean;
output.RadialOccupationError    = nr_std;

output.GibbsFit                 = GibbsFit;
output.Gibbs_z0                 = GibbsFit.z0;
output.Gibbs_z0Err              = Gz0Err;
output.Gibbs_Rs                 = GibbsFit.s*aL_um; % in um
output.Gibbs_RsErr              = GsErr*aL_um;      % in um

%% Initialize Figure

if ~isfield(opts,'Parent')
    opts.Parent = figure('color','w','Position',[300 100 1200 610],...
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

hF.UserData.InputData = input_data;
hF.UserData.OutputData = {};
    



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
%%  2D Image of Average Image
hF.UserData.Axes{1}=subplot(1,2,1,'parent',hF);
ax1= hF.UserData.Axes{1};
ax1.UserData.subplot_inds = [1 2 1];
ax1.Units='pixels';
ax1.Position=getAxesPos(1,2,1,W,H);

if size(digdata.Zdig,3)>1
    imagesc(n1,n2,imboxfilt(Zbar),'parent',ax1)
    strBoxCar = ['boxcar avg: ' num2str(opts.BinStep) '\times' num2str(opts.BinStep)];
    text(.99,.01,strBoxCar,'units','normalized','fontsize',10,...
    'verticalalignment','top','horizontalalignment','right',...
    'color','r','parent',ax1);
    ax1.XAxisLocation='Top';
    set(ax1,'FontSize',8);
else
    imagesc(n1,n2,Zbar,'parent',ax1)
end
xlabel(ax1,'position (sites)');
ylabel(ax1,'position (sites)');    
text(1,1,strRadialSummary,'units','pixels','fontsize',10,...
    'horizontalalignment','left','verticalalignment','bottom',...
    'color','r','interpreter','latex','parent',ax1);  
colormap(ax1,slanCM('magma'))
% colormap(ax1,'bone');
axis(ax1,'equal');
axis(ax1,'tight');
text(.99,.99,[strC newline strS],'units','normalized','fontsize',12,...
    'verticalalignment','top','horizontalalignment','right',...
    'color','r','parent',ax1,'interpreter','latex','fontname','arial');
ax1.XAxisLocation='Top';
ax1.YAxisLocation='Right';

set(ax1,'FontSize',8);
cc1=colorbar(ax1,'location','westoutside');
% cc1.Label.String = 'charge occupation';
set(ax1,'ydir','normal','box','on','linewidth',1);
xlim(ax1,n1_lim);
ylim(ax1,n2_lim);

%% Radial Probability Density Function (PDF)
% Average charge density

% Initialize Axis Object
hF.UserData.Axes{2}=subplot(2,2,2,'parent',hF);
ax2=hF.UserData.Axes{2};
ax2.Units='pixels';
ax2.UserData.subplot_inds = [2 2 2];
subplot_inds=ax2.UserData.subplot_inds;
ax2.Position=getAxesPos(subplot_inds(1),subplot_inds(2),subplot_inds(3),W,H);

% Plot the Data
ps(1)=histogram('BinEdges',redges,'BinCounts',nr_mean);
hold(ax2,'on')
errorbar(rcen,nr_mean,nr_std,'.','linewidth',1,'markerfacecolor',[.5 .5 .5],...
    'parent',ax2,'linestyle','none','color','k');
data_str = ['Data : $N=' num2str(round(mean(digdata.Natoms(:,1)))) '\pm' ...        
    num2str(round(std(digdata.Natoms(:,1),1))) '$~(' num2str(size(Z,3)) ' images)'];
legStr ={data_str};

% Axes Labeling
xlabel(ax2,'radial distance (sites)');
ylabel(ax2,'$n(r)$','interpreter','latex','fontsize',14);
xlim(ax2,[0 5+ceil(prctile(r_all,98))])
text(.01,1,strRadialBin,'units','normalized','horizontalalignment','left',...
    'verticalalignment','bottom','fontsize',8,'parent',ax2);
set(ax2,'box','on','linewidth',1,'fontsize',12,...
    'yaxislocation','right')

% Plot 2D Gauss Fits
rVec=linspace(0,100,100);
for nn=1:length(opts.GaussFitDensityMax)
    s       = gauss_sigma(nn,1);
    s_err   = gauss_sigma(nn,2);
    N = Iouter(nn)/nr_partial_cdf(s,rMin(nn));
    ps(end+1)=plot(rVec,nr_fit(s,rVec,N),'-');
    legStr{end+1}=gauss_str{nn};
end

% Plot Gibbs fit
ps(end+1) = plot(rVec,feval(GibbsFit,rVec),'-');
legStr{end+1} = ['Gibbs $z_0 = ' num2str(GibbsFit.z0,'%.2f') '$, $T= ' num2str(T_HOt_g,'%.1f') 't$ ' ...
    '(' num2str(T_HO_g_nK,'%.0f') ' nK)' ];
if ~isempty(opts.GaussFitDensityMax)
legend(ps,legStr,'interpreter','latex','fontsize',8,...
    'location','northeast','parent',hF);
end


    % 
    % strRadialBin = ['radial bin ' char(916) 'r:' num2str(opts.BinStep)];
    % text(.01,1,strRadialBin,'units','normalized','horizontalalignment','left',...
    %     'verticalalignment','bottom','fontsize',8,'parent',ax2);
    % 
    % rVec=linspace(0,100,100);
    % for nn=1:length(opts.GaussFitDensityMax)
    %     s=sigma_r_gauss_fit(nn);
    %     s_err=sigma_r_gauss_fit_err(nn);
    % 
    %     N = Iouter(nn)/nr_partial_cdf(s,rMin(nn));
    %     pGaussFits(nn)=plot(rVec,nr_fit(s,rVec,N),'-','Parent',ax2);
    %     legStr{nn}=['$' num2str(s*aL_um,'%.1f')  ...
    %         '(' num2str(round(10*s_err*aL_um)) ')~\mu\mathrm{m}$' ...
    %         ' $n<' num2str(opts.GaussFitDensityMax(nn)) '$'];
    % end
    % 
    % % Plot gibbs fit
    % pGaussFits(nn+1) = plot(rVec,feval(GibbsFit,rVec),'-','Parent',ax2);
    % legStr{nn+1} = ['Gibbs $z_0 = ' num2str(GibbsFit.z0,'%.1f') '$, $T= ' num2str(T_HOt_g,'%.1f') 't$ ' ...
    %     '(' num2str(T_HO_g_nK,'%.0f') ' nK)' ];
    % if ~isempty(opts.GaussFitDensityMax)
    % % legend(pGaussFits,legStr,'interpreter','latex','fontsize',8,...
    % %     'location','southeast','parent',hF);
    % legend(ax2,pGaussFits,legStr,'interpreter','latex','fontsize',8,...
    %     'location','southeast','parent',hF);
    % end
    
%% Radial Cummulative Density Function (CDF)

hpSummary = uipanel('parent',hF,'units','pixels',...
    'backgroundcolor','w','bordertype','none');
hpSummary.UserData.subplot_inds=[2 2 4];
subplot_inds=hpSummary.UserData.subplot_inds;
hpSummary.Position=getAxesPos(subplot_inds(1),subplot_inds(2),subplot_inds(3),W,H);
hF.UserData.Axes{3}=hpSummary;

ht_input = uitable('parent',hpSummary);
set(ht_input,'RowName',{},'ColumnName',{},'ColumnFormat',{'char', 'numeric'},...
    'ColumnEditable',[false true],'units','normalized',...
    'ColumnWidth',{90, 100});
ht_input.Data=input_data;
ht_input.Position(3:4) = ht_input.Extent(3:4);
ht_input.Position(1:2) = [0 1-ht_input.Extent(4)];


ht_output_1 = uitable('parent',hpSummary);
set(ht_output_1,'RowName',{},'ColumnName',{},'ColumnFormat',{'char', 'numeric'},...
    'ColumnEditable',[false false],'units','normalized',...
    'ColumnWidth',{95, 100});
ht_output_1.Data{1,1} =  ['T Gauss Equiparition [nK,t]'];
ht_output_1.Data = output_data;
ht_output_1.Position(3:4) = ht_output_1.Extent(3:4);
ht_output_1.Position(1:2) = [0 ht_input.Position(4)];

% ht_output_2 = uitable('parent',hpSummary);
% set(ht_output_2,'RowName',{},'ColumnName',{},'ColumnFormat',{'char', 'numeric'},...
%     'ColumnEditable',[false false],'units','normalized',...
%     'ColumnWidth',{90, 100});
% ht_output_2.Data{1,1} =  ['T Gauss Equiparition [nK,t]'];
% ht_output_2.Position(3:4) = ht_output_2.Extent(3:4);
% ht_output_2.Position(1:2) = [0 ht_input.Position(4)];

% hF.UserData.Axes{3}=subplot(2,2,4,'parent',hF);
% ax3=hF.UserData.Axes{3};
% ax3.Units='pixels';
% ax3.UserData.subplot_inds = [2 2 4];
% [x0, y0, w, h]=getAxesPos(2,2,4,W,H);
% ax3.Position = [x0 y0 w h];
% 
% plot(x,100*f,'k-','linewidth',2,'parent',ax3);
% xlabel(ax3,'radial distance (\mum)')
% xlim(ax3,[0 max(x)]);
% % ylim(ax3,[0 102]);
% set(ax3,'box','on','linewidth',1,'fontsize',12,...
%     'yaxislocation','right','fontname','times')
% ylabel(ax3,'% within radius')
% hold(ax3,'on');
% 
% 
% 
% plot([1 1]*rstar,[0 1]*interp1(x(2:end),100*f(2:end),rstar),'k--','linewidth',1,'parent',ax3);
% plot([rstar x(end)],[1 1]*interp1(x(2:end),100*f(2:end),rstar),'k--','linewidth',1,'parent',ax3);
% 
% if ~isempty(trap_str)    
%     text(.98,.02,trap_str,'interpreter','latex','horizontalalignment','right',...
%     'verticalalignment','bottom','fontsize',10,'units','normalized',...
%     'backgroundcolor','w','margin',1,'parent',ax3)
% end

%%
hF.UserData.Axes{1}.UserData.subplot_inds = [1 2 1];
hF.UserData.Axes{2}.UserData.subplot_inds = [2 2 2];
hF.UserData.Axes{3}.UserData.subplot_inds = [2 2 4];

hF.SizeChangedFcn=@figResize;
figResize(hF);

function figResize(src,evt)
     ax1.Position = getAxesPos(1,2,1,src.Position(3),src.Position(4));
    ax2.Position = getAxesPos(2,2,2,src.Position(3),src.Position(4));
    hpSummary.Position = getAxesPos(2,2,4,src.Position(3),src.Position(4));
ht_input.Position(3:4) = ht_input.Extent(3:4);
ht_input.Position(1:2) = [0 1-ht_input.Extent(4)];
ht_output_1.Position(3:4) = ht_output_1.Extent(3:4);
ht_output_1.Position(1:2) = [0 ht_input.Position(2)-ht_output_1.Position(4)];

end
end

function [sigmaR,sigmaR_err] = pdf_gauss2D(r)
    % Probability density function
    pdf_r = @(r,sigma) (2*pi*r.*exp(-r.^2./(2*sigma.^2))./(2*pi*sigma.^2)).*[r>0];
    % Cummulate density function (goes from 0 to 1);
    cdf_r = @(r,sigma) 1-exp(-r.^2/(2*sigma^2));
    [pdf_r_vals,pdf_r_cints] = mle(r,...
        'pdf',pdf_r,'cdf',cdf_r,...
        'Start',mean(r,'all')*sqrt(2/pi), ...
        'LowerBound',0);
    sigmaR = pdf_r_vals(1);
    sigmaR_err =(pdf_r_cints(2,1)-pdf_r_cints(1,1))*0.5;
end

function [sigmaR,sigmaR_err]=pdf_gauss_fit(r,rMin)

    % Probability density function
    pdf_r = @(r,sigma) (2*pi*r.*exp(-r.^2./(2*sigma.^2))./(2*pi*sigma.^2)).*[r>0];
    % Cummulate density function (goes from 0 to 1);
    cdf_r = @(r,sigma) (1-exp(-r.^2/(2*sigma^2))).*[r>0];


    r_expect = mean(r,'all');   % expected r
    r(r<rMin)=[];                   % remove data points smaller rMin

    % Fit it THIS ONLY WORKS WITH VERSION OF MATLAB HIGHER THAN r2019b
    try
    [pdf_r_vals,pdf_r_cints] = mle(r,...
        'pdf',pdf_r,'cdf',cdf_r,...
        'Start',r_expect*sqrt(2/pi), ...
        'LowerBound',0,'TruncationBounds',[rMin inf]);

    catch ME  
        [pdf_r_vals,pdf_r_cints] = mle(r,...
        'pdf',pdf_r,'cdf',cdf_r,...
        'Start',r_expect*sqrt(2/pi), ...
        'LowerBound',0);
    end    
    
    sigmaR = pdf_r_vals(1);
    sigmaR_err =(pdf_r_cints(2,1)-pdf_r_cints(1,1))*0.5;
end


function P=getAxesPos(num_rows,num_cols,index,figW,figH)
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

    P = [axX,axY,axWidth,axHeight];

end

