function [hF,outdata]=ixon_showBoxNumber(ixondata,xVar,opts)

global ixon_imgdir;

if nargin==2
    opts=struct;
    opts.NumberExpFit = 0;
    opts.NumberLorentzianFit=0;
    opts.NumberScale='Linear';

end


%% Sort the data by the parameter given
params=[ixondata.Params];
xvals=[params.(xVar)];

[xvals,inds]=sort(xvals,'ascend');
ixondata=ixondata(inds);
if isequal(xVar,'ExecutionDate')
    xvals=xvals-min(xvals);
end

%% Grab the box count
for kk=1:length(ixondata)
   for nn=1:size(ixondata(kk).ROI,1)
        BC=ixondata(kk).BoxCount(nn);         % Grab the box count
        Xc(kk,nn)=BC.Xc;Yc(kk,nn)=BC.Yc;        % X and Y center
        Xs(kk,nn)=BC.Xs;Ys(kk,nn)=BC.Ys;        % X and Y sigma   
        Zs(kk,nn)=BC.Ys;                          % ASSUME sZ=sY;                
        nbg(kk,nn)=BC.Nbkgd;                        % Background
        N(kk,nn)=BC.Ncounts;
   end        
end

% Convert sizes in meters
Xs = Xs*(ixondata(1).CameraInformation.PixelSize(1)*1E-6)*...
    ixondata(1).Magnification(1);
Ys = Ys*(ixondata(1).CameraInformation.PixelSize(1)*1E-6)*...
    ixondata(1).Magnification(2);


%% Outdata

outdata=struct;
outdata.xVar=xVar;
outdata.X=xvals;
outdata.N=N;

%% Exponential Decay Fit

doExpFit=opts.NumberExpFit;

if doExpFit
    myfit=fittype('A*exp(-t/tau)+B','coefficients',{'A','tau','B'},...
    'independent','t');
    opt=fitoptions(myfit);
    
    % Get some initial guesses
    tau0=max(xvals)/4;   
    tau0=.25;
    
    fout_exp={};
    for nn=1:size(N,2)  
        A0=max(N(:,nn));
        
        % Assign start point
        opt.StartPoint=[A0 tau0 min(N)];
        fout_exp{nn}=fit(xvals',N,myfit,opt);
    end
end

%% Sum E

if opts.NumberExp2SumFit
    B = 3.19e6;
        myfit=fittype(@(A1,tau1,A2,tau2,t) A1*exp(-t/tau1)+A2*exp(-t/tau2)+B,'coefficients',{'A1','tau1','A2','tau2'},...
    'independent','t');
    opt=fitoptions(myfit);
    
    % Get some initial guesses
    tau1=0.29;
    A1 = range(N);
    A2 = 5e6;
    tau2 = 100;
    
    fout_exp={};
      
    opt.StartPoint=[A1 tau1 A2 tau2];
    fout_exp{1}=fit(xvals',N,myfit,opt);
    
end

%% Lorentzian Fit

doLorentzFit=opts.NumberLorentzianFit;

if doLorentzFit && length(ixondata)>4
        X=xvals';
        Y=N;
    
    
        % Symmetric Lorentzian Decay
        myfit=fittype('A*(G/2).^2*((x-x0).^2+(G/2).^2).^(-1)+bg','coefficients',{'A','G','x0','bg'},...
            'independent','x');
        opt_fit=fitoptions(myfit);
        G0=10;
        bg=max(Y);
        A0=range(Y);
        
        inds=[Y<(bg-A0/2)];      
        
        x0=mean(X(inds));     
        opt_fit.StartPoint=[-A0 G0 x0 bg];   
        opt_fit.Robust='bisquare';
        fout_lorentz=fit(X,Y,myfit,opt_fit);
        
%         Y=N;
%         X=xvals';
% 
%         Assymmetric
%         g=@(x,a,x0,G) 2*G./(1+exp(a*(x-x0)));
%         y=@(x,a,x0,G,A,bg) A./(4*(x-x0).^2./g(x,a,x0,G).^2+1)+bg;        
%         myfit=fittype(@(a,x0,G,A,bg,x) y(x,a,x0,G,A,bg),'coefficients',{'a','x0','G','A','bg'},...
%             'independent','x'); 
%         opt_fit=fitoptions(myfit);
%         G0=10;
%         bg=min(Y);
%         A0=(max(Y)-min(Y))*(G0/2).^2;
%         inds=[Y<.8*max(Y)];
%         
%         
%         
%         x0=mean(X(inds));     
%         opt_fit.Robust='bisquare';     
% 
%         
%         opt_fit.StartPoint=[.2 x0 G0 -A0 max(Y)];  
%         fout_lorentz=fit(X,Y,myfit,opt_fit);

% 
%         XF=linspace(min(X),max(X),1000);
%         pExp=plot(XF,feval(fout_lorentz,XF),'r-','linewidth',2);
% 
%         str=['$f_0 = ' num2str(round(fout_lorentz.x0,1)) '$ kHz' newline ...
%             '$\mathrm{FWHM} = ' num2str(round(abs(fout_lorentz.G),2)) ' $ kHz'];
%         legend(pExp,{str},'interpreter','latex','location','best');
%         
% %         xlim([min(xx*1E-3) max(xx*1E-3)]);

        
    
end
%% Make Figure


% Create the name of the figure
[filepath,name,~]=fileparts(ixon_imgdir);

figDir=fullfile(ixon_imgdir,'figures');
if ~exist(figDir,'dir')
   mkdir(figDir); 
end

strs=strsplit(ixon_imgdir,filesep);
str=[strs{end-1} filesep strs{end}];

hF=figure('Name',[pad('Ixon Box Number',20) str],...
    'units','pixels','color','w','Menubar','none','Resize','off',...
    'numbertitle','off');
hF.Position(1)=1;
hF.Position(2)=50;
hF.Position(3)=500;
hF.Position(4)=400;
clf
drawnow;
uicontrol('style','text','string','iXon, box','units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',12,'fontweight','bold',...
    'position',[2 2 80 20]);

% Make axis
hax=axes;
set(hax,'box','on','linewidth',1,'fontsize',12,'units','pixels');
hold on
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
ylabel('box counts');
hax.Position(4)=hax.Position(4)-20;

for nn=1:size(ixondata(1).ROI,1)
    [cface,cedge] = ixoncolororder(nn);
   plot(xvals,N(:,nn),'o','color',cedge,'linewidth',1,'markersize',8,...
       'markerfacecolor',cface,'markeredgecolor',cedge);
end

if doExpFit
    fstrs={};
    xx=linspace(0,max(xvals),10000);
    
    for nn=1:size(N,2)
        pExp(nn)=plot(xx,feval(fout_exp{nn},xx),'r-','linewidth',1);    
        str1=['$N_0 = ' num2str(fout_exp{nn}.A,'%.2e') '$' newline ...
            '$\tau = ' num2str(round(fout_exp{nn}.tau,2)) '$' newline...
            'bg $= ' num2str(fout_exp{nn}.B,'%.2e') '$'];
        fstrs{nn}=str1;
    end   
    
    legend(pExp,fstrs,'interpreter','latex','location','best');
    hax.YLim(1)=0;
end

if opts.NumberExp2SumFit
        fstrs={};
    xx=linspace(0,max(xvals),100);
    
    for nn=1:size(N,2)
        pExp(nn)=plot(xx,feval(fout_exp{nn},xx),'r-','linewidth',1);    
%         str1=['$N_0 = ' num2str(fout_exp{nn}.A1,'%.2e') '$' newline ...
%             '$\tau = ' num2str(round(fout_exp{nn}.tau,2)) '$' newline...
%             'bg $= ' num2str(fout_exp{nn}.B,'%.2e') '$'];
        fstrs{nn}='h';
    end   
end

if doLorentzFit
    fstrs={};
    xx=linspace(0,max(xvals),100);
    
    pL=plot(xx,feval(fout_lorentz,xx),'r-','linewidth',1);    
    str1=['$x_0 = ' num2str(round(fout_lorentz.x0,1)) '$' newline ...
        '$\Gamma = ' num2str(round(fout_lorentz.G,2)) ' $'];
    fstrs={str1};
    
    legend(pL,fstrs,'interpreter','latex','location','best');
    hax.YLim(1)=0;
end

set(gca,'YScale',opts.NumberScale);


% Image directory folder string
t=uicontrol('style','text','string',str,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',6);
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];

end

