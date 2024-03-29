function [hF,outdata]=ixon_showGaussNumber(ixondata,xVar,opts)

global ixon_imgdir;

if nargin==2
    opts=struct;
    opts.NumberExpFit = 0;
    opts.NumberLorentzianFit=0;
    opts.NumberScale= 'log';
end


%% Sort the data by the parameter given
params=[ixondata.Params];
xvals=[params.(xVar)];

[xvals,inds]=sort(xvals,'ascend');
ixondata=ixondata(inds);

if isequal(xVar,'ExecutionDate')
    xvals=xvals-min(xvals);
end
%% Grab the gaussian fit outputs
for kk=1:length(ixondata)
   for nn=1:length(ixondata(kk).GaussFit)
        fout=ixondata(kk).GaussFit{nn};             % Grab the fit
        Xc(kk,nn)=fout.Xc;Yc(kk,nn)=fout.Yc;        % X and Y center
        Xs(kk,nn)=fout.Xs;Ys(kk,nn)=fout.Ys;        % X and Y sigma   
        A(kk,nn)=fout.A;                            % Amplitude
        nbg(kk,nn)=fout.nbg;                        % Background
        N(kk,nn)=2*pi*Xs(kk,nn)*Ys(kk,nn)*A(kk,nn); % Number of counts
        
        if ismember('theta',coeffnames(fout))   
            theta(kk,nn)=fout.theta;        
        end        
   end
end


% Convert sizes in meters
X_scale = (ixondata(1).CameraInformation.PixelSize(1))/...
    ixondata(1).Magnification(1);
Y_scale = (ixondata(1).CameraInformation.PixelSize(1))/...
    ixondata(1).Magnification(2);

% s1=X_scale*s1;
% s2=Y_scale*s2;



%% Outdata

outdata=struct;
outdata.xVar=xVar;
outdata.X=xvals;
outdata.N=N;

%% Exponential Decay Fit

doExpFit=opts.NumberExpFit;

if doExpFit
    myfit=fittype('A*exp(-t/tau)','coefficients',{'A','tau'},...
    'independent','t');

    myfit=fittype('A*exp(-t/tau)+B','coefficients',{'A','tau','B'},...
    'independent','t');


    opt=fitoptions(myfit);
    
    % Get some initial guesses
    tau0=max(xvals)/2;   
    tau0=.001;
    
    fout_exp={};
    for nn=1:size(N,2)  
        A0=max(N(:,nn));
        
        % Assign start point
        opt.StartPoint=[A0 tau0 0];
        fout_exp{nn}=fit(xvals',N,myfit,opt);
    end
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

hF=figure('Name',[pad('Ixon Gauss Number',20) str],...
    'units','pixels','color','w','Menubar','none','Resize','off',...
    'numbertitle','off');
hF.Position(1)=1;
hF.Position(2)=50;
hF.Position(3)=500;
hF.Position(4)=400;
clf
drawnow;
uicontrol('style','text','string','iXon, gauss','units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',12,'fontweight','bold',...
    'position',[2 2 100 20]);

% Make axis
hax=axes;
set(hax,'box','on','linewidth',1,'fontsize',12,'units','pixels');
hold on
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
ylabel('gauss counts');
hax.Position(4)=hax.Position(4)-20;

for nn=1:size(ixondata(1).ROI,1)
    [cface,cedge] = ixoncolororder(nn);
   plot(xvals,N(:,nn),'o','color',cedge,'linewidth',1,'markersize',8,...
       'markerfacecolor',cface,'markeredgecolor',cedge);
end
hax.YLim(1)=0;

if doExpFit
    fstrs={};
    xx=linspace(0,max(xvals),100);
    
    for nn=1:size(N,2)
        pExp(nn)=plot(xx,feval(fout_exp{nn},xx),'r-','linewidth',1);    
        str1=['$N_0 = ' num2str(fout_exp{nn}.A,'%.2e') '$' newline ...
            '$\tau = ' num2str(round(fout_exp{nn}.tau,4)) ' $' newline ...
            '$N_bg = ' num2str(fout_exp{nn}.B,'%.2e') '$'];
        fstrs{nn}=str1;
    end   
    
    legend(pExp,fstrs,'interpreter','latex','location','best');
    hax.YLim(1)=0;
end

if doLorentzFit
    fstrs={};
    xx=linspace(0,max(xvals),100);
    
    pL=plot(xx,feval(fout_lorentz,xx),'r-','linewidth',1);    
    str1=['$x_0 = ' num2str(round(fout_lorentz.x0,1)) '$' newline ...
        '$\Gamma = ' num2str(round(fout_lorentz.G,2)) ' $'];
    fstrs={str1};
    
    legend(pL,fstrs,'interpreter','latex','location','best');
%     hax.YLim(1)=0;
end


set(gca,'YScale',opts.NumberScale);
% Image directory folder string
t=uicontrol('style','text','string',str,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',6);
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];

end

