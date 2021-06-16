function analyzeStripes(ixondata,xVar,opts)

%% Sort the data by the parameter given
params=[ixondata.Params];
xvals=[params.(xVar)];


[xvals,inds]=sort(xvals,'ascend');
ixondata=ixondata(inds);

if nargin==2
   opts.theta=54;
   opts.rotrange=[300 350];
   opts.FitType='SmoothSquare';
%    opts.FitType='Sine';
   opts.LowThreshold=0.2;
end

%% Define the fit function

switch opts.FitType
    case 'Sine'
        foo=@(A,B,L,phi,xc,s,x) A*0.5*(1+B*sin(2*pi/L*x+phi)).*exp(-(x-xc).^2/(2*s^2));
        myfit=fittype(@(A,B,L,phi,xc,s,x) foo(A,B,L,phi,xc,s,x),...
            'independent','x','coefficients',{'A','B','L','phi','xc','s'});
        opt=fitoptions(myfit);

    case 'SmoothSquare'
        smoothSquare=@(t,delta) atan(sin(t)/delta)./atan(1/delta);
        goo=@(A,B,L,phi,xc,s,delta,x) A*0.5*(1+B*...
            smoothSquare(2*pi/L*x+phi,delta)).*...
            exp(-(x-xc).^2/(2*s^2));
        myfit=fittype(@(A,B,L,phi,xc,s,delta,x) goo(A,B,L,phi,xc,s,delta,x),...
            'independent','x','coefficients',{'A','B','L','phi','xc','s','delta'});
        opt=fitoptions(myfit);
end

% Fit tolereances
opt.Robust='bisquare';
opt.MaxIter=1000;
opt.MaxFunEvals=1000;
opt.MaxFunEvals=1000;
opt.TolFun=1E-7;

figure
tt=linspace(0,4*pi,1000);
plot(tt,smoothSquare(tt,.5));

%% Plot the figure


hF_live=figure;
hF_live.Position=[100 100 900 300];
hF_live.Color='w';

clf
ax1=subplot(121);
p1=plot(0,0,'-');
hold on
p2=plot(0,0,'r-');

xlabel('rotated vector');
ylabel('summed counts');

tt=text(.02,.98,'','units','normalized','fontsize',10,...
    'interpreter','latex','verticalalignment','top');

ax2=subplot(122);
h1=imagesc(zeros(1000,1000));
colorbar
colormap(purplemap);
caxis([0 300]);    
fRs={};
axis equal tight

Nsum=zeros(length(ixondata),1);
for kk=1:length(ixondata)
    % Grab the data and rotate
   Z=ixondata(kk).Z;
   Zrot=imrotate(Z,opts.theta,'bilinear','crop'); 
   
   % Select a sub ROI   
   Zx=sum(Zrot(200:300,:),1);   
%    Zx=sum(Zrot(:,:),1);      
   
   Nsum(kk)=sum(sum(Zrot));
   x=1:length(Zx);   
   

   xCguess=sum(Zx.*x)/sum(Zx);
   theta1=-4*pi;
   theta2=4*pi;
   Z0=max(Zx);
   L0=180;
   
   
   % Make an initial Guess
   switch opts.FitType
       case 'SmoothSquare'     
           opt.Start=[1.0*Z0 0.5 L0     0       xCguess    50      .1];     
           opt.Lower=[0.5*Z0 0.0 0      theta1  min(x)     2 0 ]; 
           opt.Upper=[1.5*Z0 1.5 L0+100 theta2  max(x)     1000 1];  
       case 'Sine'
          opt.Start=[max(Zx)     1 177 pi xCguess 70];     
           opt.Lower=[0.5*max(Zx) .5 0 -pi 100 10]; 
           opt.Upper=[max(Zx)*1.5 1.1 200 3*pi 600 200];  
   end

   % Find data point to exclue
    ilow=(Zx<Zx/max(Zx*opts.LowThreshold));   
    opt.Exclude=ilow;
   
    % Fit the data
   fRs{kk}=fit(x',Zx',myfit,opt);   
  
   % Extract parameters
  phis(kk)=mod(fRs{kk}.phi,2*pi);
   Ls(kk)=fRs{kk}.L;
   Bs(kk)=fRs{kk}.B;
   
      str=['$\phi=' num2str(phis(kk)/pi,2) '\pi$' newline ...
       '$\lambda=' num2str(Ls(kk),0) '~\mathrm{px}$' newline ...
       '$\mathrm{depth}=' num2str(round(100*Bs(kk),1)) '\%$'];
   tt.String=str;
   
   % Update the live plot
    set(p1,'XData',x,'YData',Zx);
    set(p2,'XData',x,'YData',feval(fRs{kk},x));
    set(ax1,'XLim',[1 length(x)],'YLim',[0 max(Zx)+100]);
    set(h1,'XData',1:size(Zrot,2),'YData',1:size(Zrot,1),'CData',Zrot);   
    set(ax2,'XLim',[1 size(Zrot,2)],'YLim',[1 size(Zrot,1)]);
    drawnow;
       

   
   waitforbuttonpress

   
       
end


%% Remove "bad"

i1=logical([Bs>0.3] .* [Nsum>6.5E6]');
i2=~i1;
% keyboard
%% Plot the analysis
hf2=figure;
hf2.Color='w';
hf2.Position=[100 200 600 400];
clf


subplot(221);
plot(xvals(i1),phis(i1)/pi,'o','linewidth',2)
hold on
plot(xvals(i2),phis(i2)/pi,'o','linewidth',2)
ylabel('\phi (\pi)');
xlabel(xVar);

subplot(222)
plot(xvals(i1),Ls(i1),'o','linewidth',2)
hold on
plot(xvals(i2),Ls(i2),'o','linewidth',2)
ylabel('\lambda (px)');
xlabel(xVar);

subplot(223)
plot(xvals(i1),Nsum(i1),'o','linewidth',2)
hold on
plot(xvals(i2),Nsum(i2),'o','linewidth',2)
ylabel('sum counts');
xlabel(xVar);
yL=get(gca,'YLim');
ylim([0 yL(2)]);

subplot(224)
plot(xvals(i1),Bs(i1),'o','linewidth',2)
hold on
plot(xvals(i2),Bs(i2),'o','linewidth',2)
ylabel('modulation depth');
xlabel(xVar);
% ylim([-1 1]*50+Lc);



end

