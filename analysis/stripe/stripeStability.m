function [hF,hF2] = stripeStability(stripe_data,field_gradient)

[cface1,cedge1] = ixoncolororder(1);
hF2=[];
hF=figure;
hF.Color='w';
hF.Position=[100 200 400 300];
clf

aL=532E-9;                  % Plane spacing in meters
G=field_gradient*100;   % G/cm to G/m

B0=(stripe_data.Phase/(2*pi))*G*aL;
X=stripe_data.XData;

errorbar(X,1E3*B0(:,1),B0(:,2)*1E3,'marker','o',...
    'MarkerFacecolor',cface1,'markeredgecolor',cedge1,'linestyle','none',...
    'linewidth',1.5,'color',cedge1);
hold on
ylabel('magnetic field (mG)');
xlabel(stripe_data.xVar);

if stripe_data.grabMagnetometer && isequal(stripe_data.xVar,'ExecutionDate')
    
    hF2=figure;
    hF2.Color='w';
    hF2.Position=[100 600 400 300];
    t1=stripe_data.ExecutionDate(1)/(24*60*60);
    t2=stripe_data.ExecutionDate(end)/(24*60*60);
   
   t1=datevec(t1);
   t2=datevec(t2);
   
   [T,B]=grabMagnetometer(t1,t2);  
   T=datenum(T)*24*60*60;
   T=T-min(T); 
   
   dT=T(2)-T(1);
   Nsmooth=stripe_data.Nsmooth;
   
   Tsmooth=round(dT*Nsmooth,1);
   
   B0=median(B);
   
   str=[num2str(Tsmooth) ' s smoothing time'];
   
   plot(T,smooth(B,Nsmooth),'k-','linewidth',1);
   ylabel('field (mG)');   
   xlabel('Execution Date (s)');
   
   ylim([-1 1]*5+B0);
   
   text(0.02,.98,str,'units','normalized',...
       'verticalalignment','top','fontsize',10);
end

end

