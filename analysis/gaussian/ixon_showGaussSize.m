function hF=ixon_showGaussSize(ixondata,xVar)


global ixon_imgdir


%% Sort the data by the parameter given
params=[ixondata.Params];
xvals=[params.(xVar)];

[xvals,inds]=sort(xvals,'ascend');
ixondata=ixondata(inds);

%% Grab the gaussian fit outputs
for kk=1:length(ixondata)
   for nn=1:length(ixondata(kk).GaussFit)
        fout=ixondata(kk).GaussFit{nn};             % Grab the fit
        Xc(kk,nn)=fout.Xc;Yc(kk,nn)=fout.Yc;        % X and Y center
        s1(kk,nn)=fout.s1;s2(kk,nn)=fout.s2;        % X and Y sigma   
        A(kk,nn)=fout.A;                            % Amplitude
        nbg(kk,nn)=fout.nbg;                        % Background
        N(kk,nn)=2*pi*s1(kk,nn)*s2(kk,nn)*A(kk,nn); % Number of counts
        
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

s1=X_scale*s1;
s2=Y_scale*s2;


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

hF=figure('Name',[pad('Ixon Gauss Size',20) str],...
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
xlabel(xVar,'interpreter','none','fontsize',10);
co=get(gca,'colororder');

for nn=1:size(ixondata(1).ROI,1)
    [cface,cedge] = ixoncolororder(nn);
   plot(xvals,s1(:,nn),'o-','color',cedge,'linewidth',1,'markersize',8,...
       'markerfacecolor',cface,'markeredgecolor',cedge);
end

str='$\sigma_X (\mathrm{px})$';
str='$\sigma_X (\mu\mathrm{m})$';

text(0.02,.98,str,'units','normalized','fontsize',12,'verticalalignment','cap',...
    'interpreter','latex');


% Make axis

hax2=subplot(132);
set(hax2,'box','on','linewidth',1,'fontsize',10,'units','pixels');
hold on
xlabel(xVar,'interpreter','none','fontsize',10);
co=get(gca,'colororder');
for nn=1:size(ixondata(1).ROI,1)
    [cface,cedge] = ixoncolororder(nn);
   plot(xvals,s2(:,nn),'o-','color',cedge,'linewidth',1,'markersize',8,...
       'markerfacecolor',cface,'markeredgecolor',cedge);
end
str='$\sigma_Y (\mathrm{px}^2)$';
str='$\sigma_Y (\mu\mathrm{m})$';

text(0.02,0.98,str,'units','normalized','fontsize',12,'verticalalignment','cap',...
    'interpreter','latex');



% Make axis
hax3=subplot(133);
set(hax3,'box','on','linewidth',1,'fontsize',10,'units','pixels');
hold on
xlabel(xVar,'interpreter','none','fontsize',10);
co=get(gca,'colororder');
for nn=1:size(ixondata(1).ROI,1)
    [cface,cedge] = ixoncolororder(nn);
   plot(xvals,pi*s1(:,nn).*s2(:,nn),'o-','color',cedge,'linewidth',1,'markersize',8,...
       'markerfacecolor',cface,'markeredgecolor',cedge);
end
str='$\pi \sigma_X \sigma_Y (\mathrm{px}^2)$';
str='$\pi \sigma_X \sigma_Y (\mu\mathrm{m}^2)$';
text(0.02,0.98,str,'units','normalized','fontsize',12,'verticalalignment','cap',...
    'interpreter','latex');


hax1.Position(4)=hax1.Position(4)-15;
hax2.Position(4)=hax1.Position(4);
hax3.Position(4)=hax1.Position(4);

uicontrol('style','text','string','iXon, gauss','units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',12,'fontweight','bold',...
    'position',[2 2 100 20]);

