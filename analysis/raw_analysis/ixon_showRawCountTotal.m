function hF=ixon_showRawCountTotal(ixondata,xVar,opts)
% Grab important global variables

global ixon_imgdir

params=[ixondata.Params];
acq=[ixondata.AcquisitionInformation];

xData=1:length(ixondata);
if ismember(xVar,fieldnames(params))
    xData=[params.(xVar)];
end

if ismember(xVar,fieldnames(acq))
    xData=[acq.(xVar)];
end

%% Sort the data by the parameter given

[xData,inds]=sort(xData,'ascend');
ixondata=ixondata(inds);

%% Grab data

Nsum=zeros(length(ixondata),size(ixondata(1).RawImages,3));
for kk=1:length(ixondata)
    for jj=1:size(ixondata(1).RawImages,3)
       Nsum(kk,jj)=ixondata(kk).Raw(jj).Total; 
    end
end

%% Linear Fit


if opts.FitLinear
    for jj=1:size(ixondata(1).RawImages,3)
       y=Nsum(:,jj);
       [pps{jj},pss{jj}]=polyfit(xData,y',1);
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

hF=figure('Name',['Raw Counts : '  str],...
    'units','pixels','color','w','Menubar','none','Resize','off',...
    'numbertitle','off');
hF.Position(1)=0;
hF.Position(2)=50;
hF.Position(3)=500;
hF.Position(4)=400;
clf
drawnow;



uicontrol('style','text','string','iXon','units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',12,'fontweight','bold',...
    'position',[2 2 40 20]);

% Make axis
hax=axes;
set(hax,'box','on','linewidth',1,'fontsize',14,'units','pixels');
hold on
xlabel(xVar,'interpreter','none');
ylabel('total raw counts');

hax.Position(4)=hax.Position(4)-20;

co=get(gca,'colororder');
coIN=brighten(co,.8);

legStr={};

for nn=1:size(ixondata(1).RawImages,3)
    [cface,cedge] = ixoncolororder(nn);
    if opts.FitLinear
       ps(nn)=plot(xData,polyval(pps{nn},xData),'color',cedge,'linewidth',2);
        legStr{nn}=[num2str(pps{nn}(1),'%.3e') ' counts/var'];

    end
hold on
   plot(xData,Nsum(:,nn),'o','color',cedge,'linewidth',1,'markersize',8,...
       'markerfacecolor',cface,'markeredgecolor',cedge,...
       'linewidth',2);
end
if opts.FitLinear
    legend(ps,legStr,'location','northwest');
end

% Image directory folder string
t=uicontrol('style','text','string',str,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left');
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];

end

