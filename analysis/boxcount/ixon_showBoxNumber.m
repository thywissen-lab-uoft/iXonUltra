function [hF,outdata]=ixon_showBoxNumber(ixondata,xVar,opts)

global ixon_imgdir;

if nargin==2
    opts=struct;
    opts.NumberExpFit = 0;
end


%% Sort the data by the parameter given
params=[ixondata.Params];
xvals=[params.(xVar)];

[xvals,inds]=sort(xvals,'ascend');
ixondata=ixondata(inds);

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
    myfit=fittype('A*exp(-t/tau)','coefficients',{'A','tau'},...
    'independent','t');
    opt=fitoptions(myfit);
    
    % Get some initial guesses
    tau0=max(xvals)/2;   
    
    fout_exp={};
    for nn=1:size(N,2)  
        A0=max(N(:,nn));
        
        % Assign start point
        opt.StartPoint=[A0 tau0];
        fout_exp{nn}=fit(xvals',N,myfit,opt);
    end
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
xlabel(xVar,'interpreter','none');
ylabel('box counts');
hax.Position(4)=hax.Position(4)-20;

for nn=1:size(ixondata(1).ROI,1)
    [cface,cedge] = ixoncolororder(nn);
   plot(xvals,N(:,nn),'o','color',cedge,'linewidth',1,'markersize',8,...
       'markerfacecolor',cface,'markeredgecolor',cedge);
end

if doExpFit
    strs={};
    xx=linspace(0,max(xvals),100);
    
    for nn=1:size(Natoms,2)
        pExp(nn)=plot(xx,feval(fout_exp{nn},xx),'r-','linewidth',1);    
        str=['$N_0 = ' num2str(round(fout_exp{nn}.A)) '$' newline ...
            '$\tau = ' num2str(round(fout_exp{nn}.tau,1)) ' $'];
        strs{nn}=str;
    end   
    
    legend(pExp,strs,'interpreter','latex','location','best');
    hax.YLim(1)=0;
end

% Image directory folder string
t=uicontrol('style','text','string',str,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',6);
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];

end

