function analyzeStripes(ixondata,xVar,opts)
global ixon_imgdir

%% Sort the data by the parameter given
params=[ixondata.Params];
xvals=[params.(xVar)];


[xvals,inds]=sort(xvals,'ascend');
ixondata=ixondata(inds);

if nargin==2
   opts.theta=54;
   opts.rotrange=[220 300];
   opts.FitType='Sine';
%    opts.FitType='SmoothSquare';
   opts.LowThreshold=0.2;
   opts.L0=140;
   opts.FieldGradient=210;
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

%% Prepare the live figure

hF_live=figure;
hF_live.Position=[100 100 900 300];
hF_live.Color='w';
clf

% Folder directory
strs=strsplit(ixon_imgdir,filesep);
str=[strs{end-1} filesep strs{end}];
% Image directory folder string
t=uicontrol('style','text','string',str,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',6);
t.Position(4)=t.Extent(4);
t.Position(3)=hF_live.Position(3);
t.Position(1:2)=[5 hF_live.Position(4)-t.Position(4)];

% Ixon label
uicontrol('style','text','string','iXon','units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',12,'fontweight','bold',...
    'position',[2 2 80 20]);

ax1=subplot(121);
pData=plot(0,0,'-');
hold on
pFit=plot(0,0,'r-');
% pThresh=plot(0,0,'k--');

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
   Zx=sum(Zrot(opts.rotrange(1):opts.rotrange(2),:),1);   
   
   x=1:length(Zx);      

   % Construct some initial guessing
   xCguess=sum(Zx.*x)/sum(Zx);
   theta1=-4*pi;
   theta2=4*pi;
   Z0=max(Zx);
   L0=opts.L0;   
   
   % Make an initial Guess
    switch opts.FitType
        case 'SmoothSquare'     
            opt.Start=[1.0*Z0 0.5 L0     0       xCguess    50      .1];     
            opt.Lower=[0.5*Z0 0.0 L0-50      theta1  min(x)     2 0 ]; 
            opt.Upper=[1.5*Z0 1.5 L0+50 theta2  max(x)     1000 1];  
        case 'Sine'
          opt.Start=[1.0*Z0 0.5 L0     0       xCguess    50      ];     
            opt.Lower=[0.5*Z0 0.0 L0-50      theta1  min(x)     2]; 
            opt.Upper=[1.5*Z0 1.5 L0+50 theta2  max(x)     1000];  
    end

    % Find data point to exclue
    ilow=(Zx<Zx/max(Zx*opts.LowThreshold));   
    opt.Exclude=ilow;
   
    % Fit the data
    fR=fit(x',Zx',myfit,opt);
    
    % Fit string label
    str=['$\phi=' num2str(fR.phi/pi,2) '\pi$' newline ...
       '$\lambda=' num2str(round(fR.L,1)) '~\mathrm{px}$' newline ...
       '$\mathrm{depth}=' num2str(round(100*fR.B,1)) '\%$'];   
  
   % Update the live plot
    tt.String=str;
    set(pData,'XData',x,'YData',Zx);
    set(pFit,'XData',x,'YData',feval(fR,x));
    
%     set(pThresh,'XData',[x(1) x(end)],'YData',[1 1]*opts.LowThreshold);
    
    set(ax1,'XLim',[1 length(x)],'YLim',[0 max(Zx)+100]);
    set(h1,'XData',1:size(Zrot,2),'YData',1:size(Zrot,1),'CData',Zrot);   
    set(ax2,'XLim',[1 size(Zrot,2)],'YLim',[1 size(Zrot,1)]);
    drawnow;              
    
    % Append fit object
    fRs{kk}=fR;
    Nsum(kk)=sum(sum(Zrot));
end


%% Process Fit Data
phis=zeros(length(fRs),2);
Ls=zeros(length(fRs),2);
Bs=zeros(length(fRs),2);
for kk=1:length(fRs)
    % Get Confidence interval
    ci=confint(fRs{kk});
    
    % Phase
    phis(kk,1)=mod(fRs{kk}.phi,2*pi);
    phis(kk,2)=(ci(2,4)-ci(1,4))/2;   
    
    % Wavelength
    Ls(kk,1)=fRs{kk}.L;
    Ls(kk,2)=(ci(2,3)-ci(1,3))/2;   
    
    % Modulation Depth
    Bs(kk)=fRs{kk}.B;    
    Bs(kk,2)=(ci(2,2)-ci(1,2))/2;   
end


%% Plot the analysis
[cface1,cedge1] = ixoncolororder(1);


hf2=figure;
hf2.Color='w';
hf2.Position=[100 200 600 400];
clf

% Phase
subplot(221);
errorbar(xvals,phis(:,1)/pi,phis(:,2)/pi,'marker','o',...
    'MarkerFacecolor',cface1,'markeredgecolor',cedge1,'linestyle','none',...
    'linewidth',1.5,'color',cedge1);
hold on
ylabel('phase (\pi)');
xlabel(xVar);

% Wavelength
subplot(222)
errorbar(xvals,Ls(:,1),Ls(:,2),'marker','o',...
    'MarkerFacecolor',cface1,'markeredgecolor',cedge1,'linestyle','none',...
    'linewidth',1.5,'color',cedge1);
hold on
ylabel('wavelength (px)');
xlabel(xVar);

% Total Counts
subplot(223)
plot(xvals,Nsum,'o','linewidth',2,'markerfacecolor',...
    cface1,'markeredgecolor',cedge1)
hold on
ylabel('sum counts');
xlabel(xVar);
yL=get(gca,'YLim');
ylim([0 yL(2)]);

% Modulation Depth
subplot(224)
errorbar(xvals,Bs(:,1),Bs(:,2),'marker','o',...
    'MarkerFacecolor',cface1,'markeredgecolor',cedge1,'linestyle','none',...
    'linewidth',1.5,'color',cedge1);
hold on
ylabel('modulation depth');
xlabel(xVar);
ylim([0 1]);

uicontrol('style','text','string','iXon','units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',12,'fontweight','bold',...
    'position',[2 2 80 20]);

strs=strsplit(ixon_imgdir,filesep);
str=[strs{end-1} filesep strs{end}];
% Image directory folder string
t=uicontrol('style','text','string',str,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',6);
t.Position(4)=t.Extent(4);
t.Position(3)=hf2.Position(3);
t.Position(1:2)=[5 hf2.Position(4)-t.Position(4)];

%% Phase Analysis


hf3=figure;
hf3.Color='w';
hf3.Position=[100 200 400 300];
clf

aL=532E-9;                  % Plane spacing in meters
G=opts.FieldGradient*100;   % G/cm to G/m

B0=(phis/(2*pi))*G*aL;

errorbar(xvals,1E3*B0(:,1),B0(:,2)*1E3,'marker','o',...
    'MarkerFacecolor',cface1,'markeredgecolor',cedge1,'linestyle','none',...
    'linewidth',1.5,'color',cedge1);
hold on
ylabel('magnetic field (mG)');
xlabel(xVar);


end
