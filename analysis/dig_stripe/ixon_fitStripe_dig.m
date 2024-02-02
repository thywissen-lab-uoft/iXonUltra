function [out,hF1] = ixon_fitStripe_dig(n1,n2,Zb,opts)
%ixon_fitStripe_dig Fit a binned fluorescence image to a stripe pattern.
%   When taking fluoresence images of a 2D slice of the 3D cloud,
%   application of a transverse magnetic field induces a stripe pattern to
%  form which is a sample of each plane.  This code fits this distribution
%  and also uses the fluoresence per site to describe the relative focusing
%  between each stripe.
%
%   n1   : lattice site vector which describes the columns of Zb
%   n2   : lattice site vector which describes the rows of Zb
%   Zb   : matrix of size n2xn1 which has fluoresence counts per site
%   opts : options structure
%           - SumIndex : which index are the stripes (NOT WORKING)
%           - LGuess :  wavelength guess
%           - FigNum :  figure number to assign output 
%           - ColorThreshold
if nargin~=4
    opts = struct;
end

% Default direction to look at stripes
if ~isfield(opts,'SumIndex')
    opts.SumIndex = 2;
end

% Default figure number
if ~isfield(opts,'FigNum')
    opts.FigNum = 901;
end

% Default wavelength guess is to be automatically determined
if ~isfield(opts,'LGuess')
   opts.LGuess = []; 
end

% Thresholds to throw data away and for focusing score
if ~isfield(opts,'ColorThreshold')
    opts.ColorThreshold = [1000 3000];
end

fprintf('Stripe bin fitting...')
tic
%% Fit Functions

% Transverse Fit Function
fit_exp = fittype(@(A,n0,s,n) A.*exp(-(n-n0).^2/(2*s^2)),...
    'independent','n','coefficients',{'A','n0','s'});
fit_opts_t = fitoptions(fit_exp);

% Stripe Fit Function
fit_exp_stripe = fittype(@(A,n0,s,B,L,duty,R,phi,n) ...
    A.*exp(-(n-n0).^2/(2*s^2)).*(1-B*erf_pulse_wave(L,duty,phi,R,n)),...
    'independent','n','coefficients',...
    {'A','n0','s','B','L','duty','R','phi'});
fit_opts_s = fitoptions(fit_exp_stripe);
fit_opts_s.MaxIter = 1200;
fit_opts_s.MaxFunEvals = 1e3;
fit_opts_s.Robust='bisquare';
%% Data Processing
Zb(isnan(Zb))=0;Zb(isinf(Zb))=0;        % Remove non number sites
Zb(Zb<opts.ColorThreshold(1)) = 0;           % Threshold low sites

Zs1 = sum(Zb,1);
Zs1 = Zs1(:);
n1 = n1(:);

Zs2 = sum(Zb,2);
Zs2 = Zs2(:);
n2 = n2(:);

[nn1, nn2] = meshgrid(n1,n2);
n01g = sum(Zb.*nn1,'all')/sum(Zb,'all');
n02g = sum(Zb.*nn2,'all')/sum(Zb,'all');
s1g = sqrt(sum(Zb.*(nn1-n01g).^2,'all')/sum(Zb,'all'));
s2g = sqrt(sum(Zb.*(nn2-n02g).^2,'all')/sum(Zb,'all'));

R = 3;
if opts.SumIndex == 1
    Zs = smooth(Zs1,R);
    ns = n1;
    ns0 = n01g;
    ss = s1g;

    Zt = smooth(Zs2,R);
    nt = n2;
    nt0 = n02g;
    st = s2g;
else
    Zt = smooth(Zs1,R);
    nt = n1;
    nt0 = n01g;
    st = s1g;

    Zs = smooth(Zs2,R);
    ns = n2;
    ns0 = n02g;
    ss = s2g;
end

%% Fit Transverse Distribution
fit_opts_t.StartPoint = [max(Zt) nt0 st];
fout_t = fit(nt,Zt,fit_exp,fit_opts_t);

%% Fit Stripe Distribution

% Wavelength Guess
if isempty(opts.LGuess) || isnan(opts.LGuess)
    ZsumSmooth=smooth(Zs,5);
    [yA,P]=islocalmin(ZsumSmooth,'MinSeparation',10,'MaxNumExtrema',4,...
        'MinProminence',(max(ZsumSmooth)-min(ZsumSmooth))*0.05);
    nA=diff(ns(yA));
    Pvec = movsum(P(yA),2);
    Pvec = Pvec(2:end);
    [~,ind]=max(Pvec);
    L = nA(ind);
else
    L =opts.LGuess; 
end
% Phase Guess
phiVec=linspace(0,2*pi,50);
Sphi = zeros(length(phiVec),1);
for nn=1:length(phiVec)
    phi = phiVec(nn);
    Sphi(nn) = sum(sin(pi*ns/L+phi/2).^2.*Zs,'all')/sum(Zs,'all');
end
[~,ind]=min(Sphi);
phi=phiVec(ind);

% Duty Cycle Guess
B = 0.5;

% Amplitude Guess
As = max(Zs);

% Initial Guess and bounds
fit_opts_s.StartPoint = [As ns0 ss B L 0.5 0.1 phi];
fit_opts_s.Upper = [As*1.1 ns0+10 ss*1.5 1 L*1.1 0.9 0.15 phi+pi];
fit_opts_s.Lower = [As*.8 ns0-10 ss/1.5 0.1 L/1.1 0.1 0.01 phi-pi];

% Perform the fit
fout_s = fit(ns,Zs,fit_exp_stripe,fit_opts_s);

% Stripe envelope function
stripe_envelope = @(n) fout_s.A*exp(-(n-fout_s.n0).^2/(2*fout_s.s^2));

%% Score the focusing
nL = floor(min(ns)/fout_s.L-1);
nH = ceil(max(ns)/fout_s.L+1);

seps = (nL:nH)*fout_s.L + fout_s.L*mod(fout_s.phi,2*pi)/(2*pi);
seps = round(seps);

seps(seps<min(ns))=[];
seps(seps>max(ns))=[];

[scores,centers] = ixon_stripe_dig_contrast(ns,Zb,seps,opts.ColorThreshold);
[~,ind] = max(scores);
focus_center = centers(ind);
%% Create Output Data
out = struct;
out.FitTransverse   = fout_t;
out.FitStripe       = fout_s;
out.Centers         = centers;
out.Scores          = scores;
out.Lambda          = fout_s.L;
out.Phase           = fout_s.phi;
out.Duty            = fout_s.duty;
out.ModDepth        = fout_s.B;
out.FocusCenter     = focus_center;

t = toc;
disp(['(' num2str(t,2) 's)']);
%% Plot the Results

hF1=figure(opts.FigNum);
co=get(gca,'colororder');
myc = [255,140,0]/255;
hF1.Color='w';
hF1.Position(3:4) = [770 720];
if (hF1.Position(2)+hF1.Position(4))>1000;hF1.Position(2) = 100;end
clf

% Image Plot
ax1=subplot(5,5,[1 2 3 4 6 7 8 9 11 12 13 14 16 17 18 19]);
imagesc(n1,n2,Zb);
% axis equal tight
colormap([[0 0 0];winter; [1 0 0]]);

p = get(gca,'Position');
c=colorbar('location','westoutside','fontsize',10,'color',[0 0 0 0],...
    'fontname','times');
c.Label.Color='k';
c.Label.String ='counts/site';
set(gca,'position',p)
set(gca,'ydir','normal','fontsize',10,'XAxisLocation','bottom','YColor',co(1,:),...
    'XColor',co(2,:),'fontname','times','yaxislocation','right')
caxis(opts.ColorThreshold)
hold on

% Lines for indicating the absence of atoms
for kk=1:length(seps)
    plot(get(gca,'XLim'),[1 1]*seps(kk),'--','color',myc)
end

% Text label for each string
for kk=1:length(centers)
    str = ['score=' num2str(scores(kk),'%.2e')];
    text(min(nt)+6,centers(kk),str,'horizontalalignment','left',...
        'verticalalignment','middle','fontsize',10,...
        'color',myc)
end

% Text summary of stripe fit
str = ['$\lambda=' num2str(fout_s.L,'%.2f') ',' ...
    '\phi=2\pi\cdot' num2str(round(mod(fout_s.phi,2*pi)/(2*pi),2),'%.2f') ',' ...
    '\alpha = ' num2str(round(fout_s.B,2),'%.2f') '$'];
text(5,5,str,'horizontalalignment','left',...
    'verticalalignment','bottom','fontsize',14,...
    'color',myc,'interpreter','latex','backgroundcolor',[0 0 0],'Margin',1,...
    'units','pixels')

% Plot star for most in-focused plane
plot(min(nt)+2,focus_center,'pentagram','markersize',12,...
    'markerfacecolor',myc,'markeredgecolor',myc*.8)

ax2=subplot(5,5,[5 10 15 20]);
cla
plot(Zs,ns,'k-','linewidth',1,'color','k');
hold on
nsFit = linspace(min(ns),max(ns),1e3);
plot(feval(fout_s,nsFit),nsFit,'-','color',co(1,:),'linewidth',2);
set(gca,'fontsize',12,'YAxisLocation','right','YColor',co(1,:),'fontname','times',...
    'Xaxislocation','top')
ylabel('$n_2$ (site)','interpreter','latex');
ylim([min(ns) max(ns)])
plot(stripe_envelope(nsFit),nsFit,'--','color',co(1,:),'linewidth',1)
drawnow;
for kk=1:length(seps)
    plot(get(gca,'XLim'),[1 1]*seps(kk),'--','color',myc)
end
grid on

ax3=subplot(5,5,[21 22 23 24]);
plot(nt,Zt,'k-','linewidth',1);
hold on
ntFit = linspace(min(nt),max(nt),1e3);
plot(ntFit,feval(fout_t,ntFit),'-','color',co(2,:),'linewidth',2);
set(gca,'ydir','normal','fontsize',12,'Xcolor',co(2,:),'fontname','times')
xlabel('$n_1$ (site)','interpreter','latex');
xlim([min(nt) max(nt)])
grid on

linkaxes([ax1 ax2],'y');
linkaxes([ax1 ax3],'x');

subplot(5,5,25);
plot(centers,scores,'ko','markerfacecolor',[.5 .5 .5]);
xlabel('fringe position (site)');
ylabel('score');
set(gca,'fontsize',8);
yL =get(gca,'YLim');
ylim([0 yL(2)]);
end

% rectangular pulse via erf functions
function y = erf_pulse(center,smooth_width,FWHM,xx) 
    y = 0.5.*(erf((xx+FWHM*0.5-center)./smooth_width)+erf((-xx+FWHM*0.5+center)./smooth_width));
    y = y/erf_pulse_ampl(smooth_width,FWHM);
end

% Amplitude of rectangular pulse via erf function
function y = erf_pulse_ampl(smooth_width,FWHM)
    y = 0.5.*(erf((FWHM*0.5)./smooth_width)+erf((FWHM*0.5)./smooth_width));
end

% rectangular pulse wave
function y = erf_pulse_wave(L,duty,phi,R,xx)
    x1 = min(xx);
    x2 = max(xx);
    n1 = floor(x1/L-1);
    n2 = ceil(x2/L+1);
    y = xx*0;
    for n = n1:n2 
        y = y + erf_pulse(n*L+phi/(2*pi)*L,R*L,duty*L,xx);
    end
end
