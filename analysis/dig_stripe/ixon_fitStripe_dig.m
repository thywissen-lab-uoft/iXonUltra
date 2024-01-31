function [out,hF] = ixon_fitStripe_dig(n1,n2,Zb,opts)
out = struct;
if nargin~=4
    opts = struct;
end

if ~isfield(opts,'SumIndex')
    opts.SumIndex = 2;
end

ColorThreshold = [700 3000];
%% Fit Functions

fit_exp = fittype(@(A,n0,s,n) A.*exp(-(n-n0).^2/(2*s^2)),...
    'independent','n','coefficients',{'A','n0','s'});
fit_opts_t = fitoptions(fit_exp);

fit_exp_stripe = fittype(@(A,n0,s,B,L,duty,R,phi,n) ...
    A.*exp(-(n-n0).^2/(2*s^2)).*(1-B*erf_pulse_wave(L,duty,phi,R,n)),...
    'independent','n','coefficients',{'A','n0','s','B','L','duty','R','phi'});

fit_opts_s = fitoptions(fit_exp_stripe);
%% Data Processing

Zb(isnan(Zb))=0;Zb(isinf(Zb))=0;

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

At = prctile(Zt,90);
As = prctile(Zs,90);


%% Construct Tranverse Guess

fit_opts_t.StartPoint = [At nt0 st];
fout_t = fit(nt,Zt,fit_exp,fit_opts_t);

out.FitTransverse = fout_t;

%% Construct Stripe Guess

% Wavelength Guess
ZsumSmooth=smooth(Zs,5);
% [yA,P]=islocalmax(ZsumSmooth,'MinSeparation',10,...
    % 'MaxNumExtrema',4,'MinProminence',(max(ZsumSmooth)-min(ZsumSmooth))*0.05);

[yA,P]=islocalmin(ZsumSmooth,'MinSeparation',10,...
    'MaxNumExtrema',4,'MinProminence',(max(ZsumSmooth)-min(ZsumSmooth))*0.05);

nA=diff(ns(yA));
Pvec = movsum(P(yA),2);
Pvec = Pvec(2:end);

[~,ind]=max(Pvec);

L = mean(nA);
L = nA(ind);


% Phase Guess
phiVec=linspace(0,2*pi,50);
Sphi = zeros(length(phiVec),1);
for nn=1:length(phiVec)
    phi = phiVec(nn);
    Sphi(nn) = sum(sin(pi*ns/L+phi/2).^2.*Zs,'all')/sum(Zs,'all');
end
[~,ind]=min(Sphi);
phi=phiVec(ind);

% % Duty Cycle
% dutyVec = linspace(10,90,100);
% D = zeros(length(phiVec),1);
% 
% for nn = 1:length(dutyVec)
%     D(nn) = sum(foo(ns,phi,dutyVec(nn),L,5).*Zs,'all')/sum(Zs,'all');
% end

% keyboard
B = 0.5;

fit_opts_s.StartPoint = [As ns0 ss B L phi];
fit_opts_s.StartPoint = [As ns0 ss B L 0.5 0.05 phi];

% {'A','n0','s','B','L','duty','R','phi'}

fit_opts_s.Upper = [As*1.5 ns0+10 ss*1.5 0.9 L*1.1 0.9 0.15 phi+pi];
fit_opts_s.Lower = [As/1.5 ns0-10 ss/1.5 0.1 L/1.1 0.1 0.01 phi-pi];

fit_opts_s.Robust='bisquare';
% fit_opts_s.TolFun=1e-12;
% fit_opts_s.TolX=1e-9;


fout_s = fit(ns,Zs,fit_exp_stripe,fit_opts_s)

nL = floor(min(ns)/fout_s.L-1);
nH = ceil(max(ns)/fout_s.L+1);

seps = (nL:nH)*fout_s.L + fout_s.L*mod(fout_s.phi,2*pi)/(2*pi);
seps = round(seps);

seps(seps<min(ns))=[];
seps(seps>max(ns))=[];

[scores,centers] = ixon_stripe_dig_contrast(ns,Zb,seps,ColorThreshold);

foo_env = @(n) fout_s.A*exp(-(n-fout_s.n0).^2/(2*fout_s.s^2));

out.FitStripe = fout_s;
out.Centers = centers;
out.Scores = scores;

out.Lambda = fout_s.L;
out.Phase = fout_s.phi;
out.Duty = fout_s.duty;
out.ModDepth = fout_s.B;

[~,ind] = max(scores);
out.FocusCenter = centers(ind);


%%
ca = [0 0 0];       
    cb = [0.7 .1 .6];
    cc = [linspace(ca(1),cb(1),1000)' ...
        linspace(ca(2),cb(2),1000)' linspace(ca(3),cb(3),1000)'];


hF=figure(2);
hF.Color='w';
hF.Position = [100 100 770 720];
clf
co=get(gca,'colororder');
myc = [255,140,0]/255;
subplot(5,5,[1 2 3 4 6 7 8 9 11 12 13 14 16 17 18 19]);
imagesc(n1,n2,Zb);
axis equal tight
colormap bone
% colormap(cc);

xlabel('$n_1$ (site)','interpreter','latex');
ylabel('$n_2$ (site)','interpreter','latex');
set(gca,'ydir','normal','fontsize',14,'XAxisLocation','Top','YColor',co(1,:),...
    'XColor',co(2,:),'fontname','times')
caxis(ColorThreshold)
hold on
for kk=1:length(seps)
    plot(get(gca,'XLim'),[1 1]*seps(kk),'--','color',myc)
end

for kk=1:length(centers)
    str = ['n=' num2str(centers(kk)) newline 'score=' num2str(scores(kk),'%.2e')];
    text(min(nt)+6,centers(kk),str,'horizontalalignment','left',...
        'verticalalignment','middle','fontsize',10,...
        'color',myc)
end

str = ['$\lambda=' num2str(round(fout_s.L,1)) ',' ...
    '\phi=2\pi\cdot' num2str(round(mod(fout_s.phi,2*pi)/(2*pi),2),'%.2f') ',' ...
    '\alpha = ' num2str(round(fout_s.B,2),'%.2f') '$'];

text(min(nt),min(ns),str,'horizontalalignment','left',...
    'verticalalignment','bottom','fontsize',12,...
    'color',myc,'interpreter','latex')

[bob,inds] = sort(scores,'descend');
bob = bob/max(scores);

for nn=1:length(scores)
    if bob(nn)>=0.9
        plot(min(nt)+2,centers(inds(nn)),'pentagram','markersize',12,...
            'markerfacecolor',myc,'markeredgecolor',myc*.8)
    end
end


subplot(5,5,[5 10 15 20]);
cla
plot(Zs,ns,'k-','linewidth',1,'color','k');
hold on
nsFit = linspace(min(ns),max(ns),1e3);
plot(feval(fout_s,nsFit),nsFit,'-','color',co(1,:),'linewidth',2);
set(gca,'fontsize',14,'YAxisLocation','right','YColor',co(1,:),'fontname','times')
ylabel('$n_2$ (site)','interpreter','latex');
ylim([min(ns) max(ns)])
plot(foo_env(nsFit),nsFit,'--','color',co(1,:),'linewidth',1)
drawnow;
for kk=1:length(seps)
    plot(get(gca,'XLim'),[1 1]*seps(kk),'--','color',myc)
end


subplot(5,5,[21 22 23 24]);
plot(nt,Zt,'k-','linewidth',1);
hold on
ntFit = linspace(min(nt),max(nt),1e3);
plot(ntFit,feval(fout_t,ntFit),'-','color',co(2,:),'linewidth',2);
set(gca,'ydir','normal','fontsize',14,'Xcolor',co(2,:),'fontname','times')
xlabel('$n_1$ (site)','interpreter','latex');
xlim([min(nt) max(nt)])

end

function y = erf_pulse(center,smooth_width,FWHM,xx) 
    y = 0.5.*(erf((xx+FWHM*0.5-center)./smooth_width)+erf((-xx+FWHM*0.5+center)./smooth_width));
    y = y/erf_pulse_ampl(smooth_width,FWHM);
end

function y = erf_pulse_ampl(smooth_width,FWHM)
    y = 0.5.*(erf((FWHM*0.5)./smooth_width)+erf((FWHM*0.5)./smooth_width));

end

% function y = erf_pulse_wave(smooth_width,HWHM,separation,xx)
%     x1 = min(xx);
%     x2 = max(xx);
%     n1 = floor(x1/separation-3);
%     n2 = ceil(x2/separation+3);
%     y = xx*0;
%     for n = n1:n2 
%         y = y + erf_pulse(n*separation,smooth_width,HWHM,xx);
%     end
% end

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
