function hF=ixon_fft_showBoxMoments(ixondata,xVar,opts)


global ixon_imgdir


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
        BC=ixondata(kk).fft_BoxCount(nn);         % Grab the box count
        Xc(kk,nn)=BC.Xc;Yc(kk,nn)=BC.Yc;        % X and Y center
        Xs(kk,nn)=BC.Xs;Ys(kk,nn)=BC.Ys;        % X and Y sigma   
        Zs(kk,nn)=BC.Ys;                          % ASSUME sZ=sY;                
        nbg(kk,nn)=BC.Nbkgd;                        % Background
        N(kk,nn)=BC.Ncounts;
   end        
end

%% Make Filename

% Create the name of the figure
[filepath,name,~]=fileparts(ixon_imgdir);

figDir=fullfile(ixon_imgdir,'figures');
if ~exist(figDir,'dir')
   mkdir(figDir); 
end


%% Make Figure

strs=strsplit(ixon_imgdir,filesep);
str=[strs{end-1} filesep strs{end}];

hF=figure('Name',[pad('Ixon FFT Box Moments',20) str],...
    'units','pixels','color','w','Menubar','none','Resize','off',...
    'numbertitle','off');
hF.Position(1)=500;
hF.Position(2)=50;
hF.Position(3)=800;
hF.Position(4)=400;
drawnow;

% Image directory folder string
t=uicontrol('style','text','string',str,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left');
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];

% Make axis
hax1=subplot(131);
set(hax1,'box','on','linewidth',1,'fontsize',10,'units','pixels');
hold on
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
co=get(gca,'colororder');

for nn=1:size(ixondata(1).ROI,1)
    [cface,cedge] = ixoncolororder(nn);
   plot(xvals,Xs(:,nn),'o-','color',cedge,'linewidth',1,'markersize',8,...
       'markerfacecolor',cface,'markeredgecolor',cedge);
end

str='$\sigma_X (1/\mathrm{px})$';
% str='$\sigma_X (\mu\mathrm{m})$';

text(0.02,.98,str,'units','normalized','fontsize',12,'verticalalignment','cap',...
    'interpreter','latex');


% Make axis

hax2=subplot(132);
set(hax2,'box','on','linewidth',1,'fontsize',10,'units','pixels');
hold on
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
co=get(gca,'colororder');
for nn=1:size(ixondata(1).ROI,1)
    [cface,cedge] = ixoncolororder(nn);
   plot(xvals,Ys(:,nn),'o-','color',cedge,'linewidth',1,'markersize',8,...
       'markerfacecolor',cface,'markeredgecolor',cedge);
end
str='$\sigma_Y (1/\mathrm{px})$';
% str='$\sigma_Y (\mu\mathrm{m})$';

text(0.02,0.98,str,'units','normalized','fontsize',12,'verticalalignment','cap',...
    'interpreter','latex');



% Make axis
hax3=subplot(133);
set(hax3,'box','on','linewidth',1,'fontsize',10,'units','pixels');
hold on
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
co=get(gca,'colororder');
for nn=1:size(ixondata(1).ROI,1)
    [cface,cedge] = ixoncolororder(nn);
   plot(xvals,pi*Xs(:,nn).*Ys(:,nn),'o-','color',cedge,'linewidth',1,'markersize',8,...
       'markerfacecolor',cface,'markeredgecolor',cedge);
end
str='$\pi \sigma_X \sigma_Y (1/\mathrm{px}^2)$';
% str='$\pi \sigma_X \sigma_Y (\mu\mathrm{m}^2)$';
text(0.02,0.98,str,'units','normalized','fontsize',12,'verticalalignment','cap',...
    'interpreter','latex');


hax1.Position(4)=hax1.Position(4)-15;
hax2.Position(4)=hax1.Position(4);
hax3.Position(4)=hax1.Position(4);

uicontrol('style','text','string','iXon, box','units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',12,'fontweight','bold',...
    'position',[2 2 80 20]);

