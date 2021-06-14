function analyzeStripes(ixondata,opts)

% This code analyzes the stripe patterns 

theta=53; % Rotation angle in degrees
% 'PSelect_Shim';




foo=@(A,B,L,phi,xc,s,x) A*0.5*(1+B*sin(2*pi/L*x+phi)).*exp(-(x-xc).^2/(2*s^2));
myfit=fittype(@(A,B,L,phi,xc,s,x) foo(A,B,L,phi,xc,s,x),...
    'independent','x','coefficients',{'A','B','L','phi','xc','s'});
opt=fitoptions(myfit);
opt.Robust='on';
opt.MaxIter=1000;
opt.MaxFunEvals=1000;


smoothSquare=@(t,delta) atan(sin(t)/delta)./atan(1/delta);
goo=@(A,B,L,phi,xc,s,delta,x) A*0.5*(1+B*...
    smoothSquare(2*pi/L*x+phi,delta)).*...
    exp(-(x-xc).^2/(2*s^2));

myfit2=fittype(@(A,B,L,phi,xc,s,delta,x) goo(A,B,L,phi,xc,s,delta,x),...
    'independent','x','coefficients',{'A','B','L','phi','xc','s','delta'});
opt2=fitoptions(myfit2);
opt2.Robust='bisquare';
opt2.MaxIter=1000;
opt2.MaxFunEvals=1000;
opt2.TolFun=1E-7;

opt2
% 
figure(50);
clf
ax1=subplot(121);
p1=plot(0,0,'-');
hold on
p2=plot(0,0,'r-');
xlabel('Execution Date');
ylabel('summed counts');

ax2=subplot(122);
h1=imagesc(zeros(1000,1000));
colorbar
   colormap(purplemap);
   caxis([0 300]); 

   
fRs={};

for kk=1:length(ixondata)
   Z=ixondata(kk).Z;
   Zrot=imrotate(Z,theta); 
   Zx=sum(Zrot(300:350,:),1);
%       Zx=sum(Zrot(:,:),1);
   x=1:length(Zx);   
   

   % Do the Fit
%    opt.Start=[max(Zx)     0.5 180 pi 400 100];     
%    opt.Lower=[0.5*max(Zx) 0 160 0 100 50]; 
%    opt.Upper=[max(Zx)*1.5 1 200 2*pi 600 300];  
%    fRs{kk}=fit(x',Zx',myfit,opt);        

   opt2.Start=[max(Zx)     1 177 pi 400 70 .1];     
   opt2.Lower=[0.5*max(Zx) .5 160 0 100 10 0 ]; 
   opt2.Upper=[max(Zx)*1.5 1.1 200 2*pi 600 200 1];  
   
      ilow=(Zx<Zx/max(Zx*.3));
   opt2.Exclude=ilow;
   
   fRs{kk}=fit(x',Zx',myfit2,opt2);    
   

   
   
   % Update sume plot
   set(p1,'XData',x,'YData',Zx);
    set(p2,'XData',x,'YData',feval(fRs{kk},x));

   set(ax1,'XLim',[1 length(x)],'YLim',[0 max(Zx)+100]);
   drawnow;   
   
   
  % Update image 
   set(h1,'XData',1:size(Zrot,2),'YData',1:size(Zrot,1),'CData',Zrot);   
   set(ax2,'XLIm',[1 size(Zrot,2)],'YLim',[1 size(Zrot,1)]);
    drawnow;
   
    
   phis(kk)=fRs{kk}.phi;
%       phis(kk)=fRs{kk}.dx;

   Ls(kk)=fRs{kk}.L;
%    waitforbuttonpress
end


figure(1);
clf
subplot(131);
plot(phis/pi,'o-','linewidth',2)
ylabel('\phi (\pi)');
xlabel('executio date');

subplot(132);

% Plot a circle.
tt = linspace(0, 2*pi, 100);

% Plot circle.
plot(cos(tt),sin(tt), 'k-', 'LineWidth', 2);

% Plot center.
hold on;
plot(cos(phis),sin(phis), 'ko', 'LineWidth', 2, 'MarkerSize', 8);


subplot(133)
plot(Ls,'linewidth',2)
ylabel('\lambda (px)');
xlabel('execution date');

% ylim([0 200]);





end

