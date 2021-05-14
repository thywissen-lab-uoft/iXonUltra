function hF=ixon_showGaussStats(ixondata)


global ixon_imgdir;

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

% s1=X_scale*s1;
% s2=Y_scale*s2;



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

hF=figure('Name',[pad('Gauss Stats',20) str],...
    'units','pixels','color','w','Menubar','none','Resize','off',...
    'numbertitle','off');
hF.Position(1)=50;
hF.Position(2)=50;
hF.Position(3)=600;
hF.Position(4)=400;
drawnow;

% Image directory folder string
t=uicontrol('style','text','string',str,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',6);
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];

subplot(231);
histogram(N,20)
set(gca,'box','on','linewidth',1,'fontsize',12)
% text(5,5,'number','interpreter','latex',...
%     'units','pixels','backgroundcolor',[1 1 1],'verticalalignment','bottom',...
%     'edgecolor','k','margin',1);
xlabel('counts','fontsize',8)

subplot(232);
histogram(nbg,20)
set(gca,'box','on','linewidth',1,'fontsize',12)
% text(5,5,'bkgd','interpreter','latex',...
%     'units','pixels','backgroundcolor',[1 1 1],'verticalalignment','bottom',...
%     'edgecolor','k','margin',1);
xlabel('bkgd','fontsize',8)


subplot(233);
histogram(Xc,20)
set(gca,'box','on','linewidth',1,'fontsize',12)
% text(5,5,'$X_c ~(\mathrm{px})$','interpreter','latex',...
%     'units','pixels','backgroundcolor',[1 1 1],'verticalalignment','bottom',...
%     'edgecolor','k','margin',1);
xlabel('x center (px)','fontsize',8)


subplot(234);
histogram(Yc,20)
set(gca,'box','on','linewidth',1,'fontsize',12)
% text(5,5,'$Y_c ~(\mathrm{px})$','interpreter','latex',...
%     'units','pixels','backgroundcolor',[1 1 1],'verticalalignment','bottom',...
%     'edgecolor','k','margin',1);
xlabel('y center (px)','fontsize',8)

subplot(235);
histogram(s1,20)
set(gca,'box','on','linewidth',1,'fontsize',12)
% text(5,5,'$\sigma_X ~(\mu\mathrm{m})$','interpreter','latex',...
%     'units','pixels','backgroundcolor',[1 1 1],'verticalalignment','bottom',...
%     'edgecolor','k','margin',1);
xlabel('\sigma1 (px)','fontsize',8)

subplot(236);
histogram(s2,20)
set(gca,'box','on','linewidth',1,'fontsize',12)
% text(5,5,'$\sigma_Y ~(\mu\mathrm{m})$','interpreter','latex',...
%     'units','pixels','backgroundcolor',[1 1 1],'verticalalignment','bottom',...
%     'edgecolor','k','margin',1);
xlabel('\sigma2 (px)','fontsize',8)

end

