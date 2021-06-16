function [hF,out]=ixon_showGaussAspectRatio(ixondata,xVar)
% Grab important global variables

 global ixon_imgdir;


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

%% out data

out=struct;
out.s1=s1;
out.s2=s2;
out.s1s2=s1./s2;

%% Make Figure

% Create the name of the figure
[filepath,name,~]=fileparts(ixon_imgdir);

figDir=fullfile(ixon_imgdir,'figures');
if ~exist(figDir,'dir')
   mkdir(figDir); 
end

strs=strsplit(ixon_imgdir,filesep);
str=[strs{end-1} filesep strs{end}];

hF=figure('Name',[pad('Ixon Gauss Aspect Ratio',20) str],...
    'units','pixels','color','w','Menubar','none','Resize','off',...
    'numbertitle','off');
hF.Position(1)=0;
hF.Position(2)=480;
hF.Position(3)=400;
hF.Position(4)=400;
drawnow;

% Image directory folder string
t=uicontrol('style','text','string',str,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',6);
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];

uicontrol('style','text','string','iXon, gauss','units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',12,'fontweight','bold',...
    'position',[2 2 100 20]);

% Make axis
hax=axes;
set(hax,'box','on','linewidth',1,'fontsize',12,'units','pixels');
hold on
xlabel(xVar,'interpreter','none');
ylabel('aspect ratio \sigma_1/\sigma_2');

hax.Position(4)=hax.Position(4)-20;

co=get(gca,'colororder');

for nn=1:size(ixondata(1).ROI,1)
    [cface,cedge] = ixoncolororder(nn);    
   plot(xvals,out.s1s2(:,nn),'o','color',cedge,'linewidth',1,'markersize',8,...
       'markerfacecolor',cface,'markeredgecolor',cedge);
end



end

