function BinStripe = bin_StripeFitCircular(n1,n2,Zb,opts)
%bin_StripeFit Fit a binned fluorescence image to a stripe pattern.
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
%           - Threshold
if nargin~=4
    opts = struct;
end

% Default direction to look at stripes
if ~isfield(opts,'SumIndex')
    opts.SumIndex = 1; % 1 for trees % 2 for fallen trees
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
if ~isfield(opts,'Threshold')
    opts.Threshold = [1000 3000];

end

[nn1,nn2]=meshgrid(n1,n2);


%% Circular Wave Phase
% This function creates a phase map
% L = wavelength
% R      = radius of curvature
% theta  = angle
% phi0   = phase offset

n10 = mean(n1);
n20 = mean(n2);
y0=n20;


x0=n10;
x0=0;


theta0=atan2(n20,n10);
phase_func = @(L,R,theta,phi0,xx,yy) ...
    2*pi/L*(R^2+2*R*(cos(theta)*(xx-x0)+sin(theta)*(yy-y0))+((xx-x0).^2+(yy-y0).^2)).^(1/2)-2*pi*R/L+phi0;



keyboard
% Lguess = [

%% Fit Functions

% boop=
% 
% 
% phase_func = @(n1,n2) 
% 
% 
% 
% % Transverse Fit Function
% fit_exp = fittype(@(A,n0,s,n) A.*exp(-(n-n0).^2/(2*s^2)),...
%     'independent','n','coefficients',{'A','n0','s'});
% fit_opts_t = fitoptions(fit_exp);
% 
% % Stripe Fit Function
% myfunc = @(A,n0,s,B,L,duty,R,phi,n) ...
%     A.*exp(-(n-n0).^2/(2*s^2)).*(1-B*erf_pulse_wave(L,duty,phi,R,n));
% 
% % Stripe Fit Function
% % myfunc = @(A,n0,s,B,L,duty,R,phi,n) ...
% %     A.*exp(-(n-n0).^2/(2*s^2)).*(B+(1-B)*erf_pulse_wave(L,duty,phi,R,n));
% 
% fit_exp_stripe = fittype(@(A,n0,s,B,L,duty,R,phi,n) myfunc(A,n0,s,B,L,duty,R,phi,n),...
%     'independent','n','coefficients',...
%     {'A','n0','s','B','L','duty','R','phi'});
% fit_opts_s = fitoptions(fit_exp_stripe);
% fit_opts_s.MaxIter = 1200;
% fit_opts_s.MaxFunEvals = 1e3;
% 
% 
% % Stripe Amplitude Function
% phase_func = @(L,phi,site) 2*pi*(site/L)+phi+pi/2;
% 
% 
% % fit_opts_s.Robust='bisquare';
% %% Data Processing
% Zb(isnan(Zb))=0;
% Zb(isinf(Zb))=0;        % Remove non number sites
% 
% 
% 
% Zbcopy=Zb;
% 
% Zbcopy(Zbcopy<opts.Threshold(1)) = 0;           % Threshold low sites
% 
% Zs1 = sum(Zbcopy,1);
% Zs1 = Zs1(:);
% n1 = n1(:);
% 
% Zs2 = sum(Zbcopy,2);
% Zs2 = Zs2(:);
% n2 = n2(:);
% 
% 
% [nn1, nn2] = meshgrid(n1,n2);
% n01g = sum(Zbcopy.*nn1,'all')/sum(Zbcopy,'all');
% n02g = sum(Zbcopy.*nn2,'all')/sum(Zbcopy,'all');
% s1g = sqrt(sum(Zbcopy.*(nn1-n01g).^2,'all')/sum(Zbcopy,'all'));
% s2g = sqrt(sum(Zbcopy.*(nn2-n02g).^2,'all')/sum(Zbcopy,'all'));
% 
% R = 3;
% if opts.SumIndex == 1
%     Zs = smooth(Zs1,R);
%     ns = n1;
%     ns0 = n01g;
%     ss = s1g;
%     Zt = smooth(Zs2,R);
%     nt = n2;
%     nt0 = n02g;
%     st = s2g;
% else
%     Zt = smooth(Zs1,R);
%     nt = n1;
%     nt0 = n01g;
%     st = s1g;
% 
%     Zs = smooth(Zs2,R);
%     ns = n2;
%     ns0 = n02g;
%     ss = s2g;
% end
% 
% %% Fit Transverse Distribution
% fit_opts_t.StartPoint = [max(Zt) nt0 st];
% [fout_t,gof_t] = fit(nt,Zt,fit_exp,fit_opts_t);
% %% Fit Stripe Distribution
% 
% % Wavelength Guess
% if isempty(opts.LGuess) || isnan(opts.LGuess)
%     ZsumSmooth=smooth(Zs,5);
%     [yA,P]=islocalmin(ZsumSmooth,'MinSeparation',10,'MaxNumExtrema',4,...
%         'MinProminence',(max(ZsumSmooth)-min(ZsumSmooth))*0.05);
%     nA=diff(ns(yA));
%     Pvec = movsum(P(yA),2);
%     Pvec = Pvec(2:end);
%     [~,ind]=max(Pvec);
%     L = nA(ind);
% else
%     L =opts.LGuess; 
% end
% 
% % Amplitude Guess
% As = max(Zs);
% 
% % Phase Guess
% phiVec=linspace(0,2*pi,50);
% Ephi = zeros(length(phiVec),1);
% for nn=1:length(phiVec)
%     phi = phiVec(nn);    
%     y = myfunc(As,ns0,ss,0.8,L,0.5,0.05,phi,ns);
%     Ephi(nn)=sum((y-Zs).^2);
% end
% [~,ind]=min(Ephi);
% phi=phiVec(ind);
% 
% B_list = linspace(0,1,20);
% D_list = linspace(.1,.9,20);
% [BB,DD]=meshgrid(B_list,D_list);
% BB=BB(:);
% DD=DD(:);
% EE=zeros(length(BB),1);
% for kk=1:length(BB)
%     y = myfunc(As,ns0,ss,BB(kk),L,DD(kk),1,phi,ns);
%     EE(kk)=sum((y-Zs).^2);
% end
% [~,ind] = min(EE);
% B = BB(ind);
% D = DD(ind);
% 
% % Initial Guess and bounds
% fit_opts_s.StartPoint = [As ns0 ss B L D .5 phi];
% fit_opts_s.Upper = [As*2 ns0+10 (ss*1.5+10) 1 L+2 1 5 phi+1.5*pi];
% fit_opts_s.Lower = [As*.5 ns0-10 0 0.01 L-2 0 0.1 phi-1.5*pi];
% 
% % Perform the fit 
% try
%     [fout_s ,gof_s] = fit(ns,Zs,fit_exp_stripe,fit_opts_s);
% end
% 
% % disp([fout_s.A fout_s.n0 fout_s.s fout_s.B fout_s.L fout_s.duty fout_s.R fout_s.phi] - ...
%     % fit_opts_s.StartPoint);
% % 
% % fprintf('dA');fprintf('%.0e',[fout_s.A - As]);fprintf(', ');
% % fprintf('dn');fprintf('%.1f',[fout_s.n0 - ns0]);fprintf(', ');
% % fprintf('dB');fprintf('%.2f',[fout_s.B - B]);fprintf(', ');
% % fprintf('dL');fprintf('%.1f',[fout_s.L - L]);fprintf(', ');
% % fprintf('dD');fprintf('%.1f',[fout_s.duty - D]);fprintf(', ');
% % fprintf('dR');fprintf('%.1f',[fout_s.R - 0.5]);fprintf(', ');
% % fprintf('dP');fprintf('%.1f',[fout_s.phi - phi]);fprintf(')');
% 
% % disp(fout_s. - As)
% % 
% % {'A','n0','s','B','L','duty','R','phi'}
% % Stripe envelope function
% stripe_envelope = @(n) fout_s.A*exp(-(n-fout_s.n0).^2/(2*fout_s.s^2));
% 
% %% Score the focusing
% nL = floor(min(ns)/fout_s.L-1);
% nH = ceil(max(ns)/fout_s.L+1);
% 
% seps = (nL:nH)*fout_s.L + fout_s.L*mod(fout_s.phi,2*pi)/(2*pi);
% seps = round(seps);
% 
% seps(seps<min(ns))=[];
% seps(seps>max(ns))=[];
% 
% [scores,centers] = bin_StripeScore(ns,Zbcopy,seps,opts.Threshold,opts.SumIndex);
% [~,ind] = max(scores);
% focus_center = centers(ind);
% 
% [scores_sort,inds] = sort(scores,'descend');
% scores_sort = scores_sort(1:3);
% scores_sort = scores_sort/max(scores_sort);
% centers_sort = centers(inds(1:3));
% 
% scores_sort(isnan(scores_sort))=0;
% 
% 
% 
% myfit = fittype('-A*(x-x0).^2+B','independent','x','coefficients',{'A','x0','B'});
% 
% fitopt=fitoptions(myfit);
% fitopt.Start = [1e-3 centers_sort(1) 1];
% fitopt.Upper = [1 max(centers_sort)+5 1.1];
% fitopt.Lower = [0 min(centers_sort)-5 .9];
% 
% 
%  fout_centers = fit(centers_sort,scores_sort,myfit,fitopt);
%  FocusCenterFit=fout_centers.x0;
%  %% Score Try 2
% %  
% %  if opts.SumIndex ==1     
% %      site_map = nn1;
% %  else
% %      site_map = nn2;
% %  end
% %  
% %  phase_map = phase_func(fout_s.L,fout_s.phi,site_map);
% % 
% %  fringe_map = round(phase_map/(2*pi));
% %  
% %  fringe_vals = unique(fringe_map);
% %  fringe_locs = zeros(length(fringe_vals),1);
% %  fringe_scores = zeros(length(fringe_vals),1);
% %  fringe_coms = zeros(length(fringe_vals),1);
% %  
% %  % Score Via kmeans clustering
% %  for kk=1:length(fringe_vals)
% %      this_fringe = [fringe_map==fringe_vals(kk)];     
% %      loc = sum(this_fringe.*site_map,'all')/sum(this_fringe,'all');
% %      fringe_locs(kk)=loc;        
% %     Zsub = Zb.*this_fringe;
% %     fringe_coms(kk) = sum(Zbcopy.*site_map.*this_fringe,'all')/sum(Zbcopy.*this_fringe,'all');
% % 
% %      
% %      if sum(Zsub,'all')/sum(Zb,'all')>0.15 && ~isnan(fringe_coms(kk))
% %          data = Zsub(:);
% %          data(data==0)=[];         
% %          [idx,C,sumd,D] = kmeans(data,2);
% %          fringe_scores(kk) = max(C);  
% %      else
% %          fringe_scores(kk)=0;
% %      end   
% %  end
% %  
% %   X              = fringe_locs;
% %   X              = fringe_coms;
% %  
% %  Y              = fringe_scores;
% %  
% %  bad_inds       = [fringe_scores==0];
% %  
% %  X(bad_inds)    = [];
% %  Y(bad_inds)    = [];
% %    
% %  [Y0,ind]       = max(Y); 
% %  X0             = X(ind);
% %  
% % myfit = fittype('-A*(x-x0).^2+B','independent','x','coefficients',{'A','x0','B'});
% % fitopt=fitoptions(myfit);
% % fitopt.Start = [1e-3 X0 1];
% % fitopt.Upper = [1 X0+20 1.1];
% % fitopt.Lower = [0 X0-20 .9];
% % 
% % 
% % fout_focus = fit(X,Y/Y0,myfit,fitopt);
% % if X0<90
% %     keyboard
% % end
% % focus_center = X0;
% % FocusCenterFit=fout_focus.x0;
% 
% % keyboard
% 
% %% Create Output Data
% BinStripe = struct;
% BinStripe.FitTransverse       = fout_t;
% BinStripe.FitStripe           = fout_s;
% BinStripe.RSquareStripe       = gof_s.rsquare;
% BinStripe.RSquareTransverse   = gof_t.rsquare;
% BinStripe.Lambda              = fout_s.L;
% BinStripe.Phase               = fout_s.phi;
% BinStripe.Duty                = fout_s.duty;
% BinStripe.ModDepth            = fout_s.B;
% 
% 
% BinStripe.Separations         = seps;
% BinStripe.Centers             = centers;
% BinStripe.Scores              = scores;
% 
% BinStripe.FocusCenter         = focus_center;
% BinStripe.FocusCenterFit      = FocusCenterFit;
% 
% BinStripe.SumIndex              = opts.SumIndex;
% BinStripe.Site2Phase          = @(site) mod(phase_func(fout_s.L,fout_s.phi,site)+pi,2*pi)-pi;
% 
% 
% 
% 
% BinStripe.Counts              = sum(Zbcopy,'all');
% BinStripe.Opts                = opts;
% % keyboard
% %% Plot the Results
% %{
% hF1=figure(opts.FigNum);
% co=get(gca,'colororder');
% myc = [255,140,0]/255;
% hF1.Color='w';
% hF1.Position(3:4) = [770 720];
% if (hF1.Position(2)+hF1.Position(4))>1000;hF1.Position(2) = 100;end
% clf
% 
% % Image Plot
% ax1=subplot(5,5,[1 2 3 4 6 7 8 9 11 12 13 14 16 17 18 19]);
% imagesc(n1,n2,Zb);
% % axis equal tight
% colormap([[0 0 0];winter; [1 0 0]]);
% 
% p = get(gca,'Position');
% c=colorbar('location','westoutside','fontsize',10,'color',[0 0 0 0],...
%     'fontname','times');
% c.Label.Color='k';
% c.Label.String ='counts/site';
% set(gca,'position',p)
% set(gca,'ydir','normal','fontsize',10,'XAxisLocation','bottom','YColor',co(1,:),...
%     'XColor',co(2,:),'fontname','times','yaxislocation','right')
% caxis(opts.ColorThreshold)
% hold on
% 
% % Lines for indicating the absence of atoms
% for kk=1:length(seps)
%     plot(get(gca,'XLim'),[1 1]*seps(kk),'--','color',myc)
% end
% 
% % Text label for each string
% for kk=1:length(centers)
%     str = ['score=' num2str(scores(kk),'%.2e')];
%     text(min(nt)+6,centers(kk),str,'horizontalalignment','left',...
%         'verticalalignment','middle','fontsize',10,...
%         'color',myc)
% end
% 
% % Text summary of stripe fit
% str = ['$\lambda=' num2str(fout_s.L,'%.2f') ',' ...
%     '\phi=2\pi\cdot' num2str(round(mod(fout_s.phi,2*pi)/(2*pi),2),'%.2f') ',' ...
%     '\alpha = ' num2str(round(fout_s.B,2),'%.2f') '$'];
% text(5,5,str,'horizontalalignment','left',...
%     'verticalalignment','bottom','fontsize',14,...
%     'color',myc,'interpreter','latex','backgroundcolor',[0 0 0],'Margin',1,...
%     'units','pixels')
% 
% % Plot star for most in-focused plane
% try
% plot(min(nt)+2,focus_center,'pentagram','markersize',12,...
%     'markerfacecolor',myc,'markeredgecolor',myc*.8)
% end
% 
% ax2=subplot(5,5,[5 10 15 20]);
% cla
% plot(Zs,ns,'k-','linewidth',1,'color','k');
% hold on
% nsFit = linspace(min(ns),max(ns),1e3);
% plot(feval(fout_s,nsFit),nsFit,'-','color',co(1,:),'linewidth',2);
% set(gca,'fontsize',12,'YAxisLocation','right','YColor',co(1,:),'fontname','times',...
%     'Xaxislocation','top')
% ylabel('$n_2$ (site)','interpreter','latex');
% ylim([min(ns) max(ns)])
% plot(stripe_envelope(nsFit),nsFit,'--','color',co(1,:),'linewidth',1)
% drawnow;
% for kk=1:length(seps)
%     plot(get(gca,'XLim'),[1 1]*seps(kk),'--','color',myc)
% end
% grid on
% strR = ['$R^2 = ' num2str(gof_t.rsquare,'%.3f') '$'];
% text(4,2,strR,'fontsize',10,'interpreter','latex','verticalalignment','bottom',...
%     'units','pixels');
% 
% ax3=subplot(5,5,[21 22 23 24]);
% plot(nt,Zt,'k-','linewidth',1);
% hold on
% ntFit = linspace(min(nt),max(nt),1e3);
% plot(ntFit,feval(fout_t,ntFit),'-','color',co(2,:),'linewidth',2);
% set(gca,'ydir','normal','fontsize',12,'Xcolor',co(2,:),'fontname','times')
% xlabel('$n_1$ (site)','interpreter','latex');
% xlim([min(nt) max(nt)])
% grid on
% strR = ['$R^2 = ' num2str(gof_s.rsquare,'%.3f') '$'];
% text(4,2,strR,'fontsize',10,'interpreter','latex','verticalalignment','bottom',...
%     'units','pixels')
% 
% linkaxes([ax1 ax2],'y');
% linkaxes([ax1 ax3],'x');
% 
% subplot(5,5,25);
% plot(centers,scores,'ko','markerfacecolor',[.5 .5 .5]);
% xlabel('fringe position (site)');
% ylabel('score');
% set(gca,'fontsize',8);
% yL =get(gca,'YLim');
% ylim([0 yL(2)]);
% xlim([min(n1) max(n1)])
%}
end

% % rectangular pulse via erf functions
% function y = erf_pulse(center,smooth_width,FWHM,xx) 
%     y = 0.5.*(erf((xx+FWHM*0.5-center)./smooth_width)+erf((-xx+FWHM*0.5+center)./smooth_width));
%     y = y/erf_pulse_ampl(smooth_width,FWHM);
% end
% 
% % Amplitude of rectangular pulse via erf function
% function y = erf_pulse_ampl(smooth_width,FWHM)
%     y = 0.5.*(erf((FWHM*0.5)./smooth_width)+erf((FWHM*0.5)./smooth_width));
% end
% 
% % rectangular pulse wave
% function y = erf_pulse_wave(L,duty,phi,R,xx)
%     x1 = min(xx);
%     x2 = max(xx);
%     n1 = floor(x1/L-1);
%     n2 = ceil(x2/L+1);
%     y = xx*0;
%     for n = n1:n2 
%         y = y + erf_pulse(n*L+phi/(2*pi)*L,R,duty*L,xx);
%     end
% end
