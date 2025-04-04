function [hF1,hF2,hF3] = fieldAnalysis(stripe_data,opts)

[cface1,cedge1] = ixoncolororder(1);


aL=opts.LatticeSpacing;                  % Plane spacing in meters
G=opts.FieldGradient*100;                   % G/cm to G/m

B0=(stripe_data.Phase/(2*pi))*G*aL;
B0(:,1)=B0(:,1)-B0(1,1); % Initialize at field =0

mag=83; % assumed microscope magnification (16 um/pix --> 16/83 um px
pxsize=16/83;



X=stripe_data.XData;


%% Field with modulus

hF1=figure;
hF1.Color='w';
hF1.Position=[100 200 600 300];
clf
errorbar(X,1E3*B0(:,1),B0(:,2)*1E3,'marker','o',...
    'MarkerFacecolor',cface1,'markeredgecolor',cedge1,'linestyle','none',...
    'linewidth',1.5,'color',cedge1);
hold on
ylabel('magnetic field (mG)');
xlabel(stripe_data.xVar);

%% Field without modulus
hF2=figure;
hF2.Color='w';
hF2.Position=[100 200 600 300];
clf
errorbar(X,1E3*mod(B0(:,1),G*aL),B0(:,2)*1E3,'marker','o',...
    'MarkerFacecolor',cface1,'markeredgecolor',cedge1,'linestyle','none',...
    'linewidth',1.5,'color',cedge1);
hold on
ylabel('magnetic field (mG)');
xlabel(stripe_data.xVar);
xlim([min(X) max(X)]);
ylim([-.1 1.1]*G*aL*1E3);

%% Field without modulus
hF1=figure;
hF1.Color='w';
hF1.Position=[100 200 600 300];
clf

pD=errorbar(X,1E3*B0(:,1),B0(:,2)*1E3,'marker','o',...
    'MarkerFacecolor',cface1,'markeredgecolor',cedge1,'linestyle','none',...
    'linewidth',1.5,'color',cedge1);
hold on
ylabel('magnetic field (mG)');
xlabel(stripe_data.xVar);
xlim([min(X) max(X)]);


if isequal(opts.FitType,'Exp')
    Bdata=B0(:,1);
    
    myfit=fittype('B1*exp(-t/tau)+B0','coefficients',{'B1','B0','tau'},...
        'independent','t');
    fitopt=fitoptions(myfit);
    fitopt.Start=[range(Bdata) Bdata(end) mean(X)];
    
    fR=fit(X,Bdata,myfit,fitopt);
    
    xx=linspace(min(X),max(X),1000);
    
    
    pF=plot(xx,feval(fR,xx)*1e3,'r-','linewidth',2);
    hold on    
    
    str=['$\tau=' num2str(round(fR.tau,1)) '$'];
    legend(pF,{str},'interpreter','latex','location','best');
end



%%



if stripe_data.grabMagnetometer && isequal(stripe_data.xVar,'ExecutionDate')
   
    hF3=figure;
    hF3.Color='w';
    hF3.Position=[100 600 600 300];
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
   text(0.02,.98,str,'units','normalized',...
       'verticalalignment','top','fontsize',10);
   xlim([min(X) max(X)]);

else
    hF3=[];
end

end

function [fitFunc,outparams]=modExpFit(x,y,Bp,pGuess)


fitFunc=@ (p,t) mod(p(1)*exp(-t/p(2))+p(3),Bp);



% xx=linspace(0,100,1000);

% da=atan2(sin(xx), cos(xx))

% figure
% % plot(xx/(pi),abs(da))

% err_func = @(p) norm(mod(y-myfunc(p,x),Bp));

% err_func = @(p) norm(angdiff(mod(fitFunc(p,x),Bp),y));


err_func = @(p) norm(atan2(sin(fitFunc(p,x)-y),cos(fitFunc(p,x)-y)));
% keyboard
options=optimset('MaxFunEvals', 1000000, ...
    'MaxIter',1000000, 'Display', 'on', 'TolX', 1e-12);

[x,fval,exitflag,output] = fminsearch(err_func,pGuess,options);

outparams=x;




end
