function analyzeStripes(ixondata,opts)

% This code analyzes the stripe patterns 

theta=53; % Rotation angle in degrees
% 'PSelect_Shim';

foo=@(A,B,L,phi,xc,s,x) A*0.5*(1+B*sin(2*pi/L*x+phi)).*exp(-(x-xc).^2/(2*s^2));
myfit=fittype(@(A,B,L,phi,xc,s,x) foo(A,B,L,phi,xc,s,x),...
    'independent','x','coefficients',{'A','B','L','phi','xc','s'});

% foo=@(A,B,L,dx,xc,s,x) A*0.5*(1+B*sin(2*pi/L*(x+dx-50))).*exp(-(x-xc).^2/(2*s^2));
% myfit=fittype(@(A,B,L,dx,xc,s,x) foo(A,B,L,dx,xc,s,x),...
%     'independent','x','coefficients',{'A','B','L','dx','xc','s'});


opt=fitoptions(myfit);



figure(50);
clf

fRs={};

for kk=1:length(ixondata)
   Z=ixondata(kk).Z;
   Zrot=imrotate(Z,theta);
   
   
   Zx=sum(Zrot(300:350,:),1);
%       Zx=sum(Zrot(:,:),1);

   x=1:length(Zx);   
%    
   opt.Start=[max(Zx)     0.5 180 pi 400 100];   
%    opt.Lower=[0.5*max(Zx) 0 50 0 100 50]; 
%    opt.Upper=[max(Zx)*1.5 1 250 2*pi 600 300];   
   
%     opt.Start=[max(Zx)     0.5 180 100 400 100];   
%    opt.Lower=[0.5*max(Zx) 0 50 0 100 50]; 
%    opt.Upper=[max(Zx)*1.5 1 250 200 600 300]; 
   
   
   fRs{kk}=fit(x',Zx',myfit,opt);     
   
%    disp(fR);   
   
   subplot(121);
   plot(x,Zx);
   hold on
   plot(fRs{kk});
   
   subplot(122);
   imagesc(Zrot);
   colormap(purplemap);
   caxis([0 300]);
%    waitforbuttonpress
   
   clf 
    
   phis(kk)=fRs{kk}.phi;
%       phis(kk)=fRs{kk}.dx;

   Ls(kk)=fRs{kk}.L;
   
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

