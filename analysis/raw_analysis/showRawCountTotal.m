function hF=showRawCountTotal(atomdata,xVar,opts)
% Grab important global variables

global imgdir

params=[atomdata.Params];
acq=[atomdata.AcquisitionInformation];

xData=1:length(atomdata);
if ismember(xVar,fieldnames(params))
    xData=[params.(xVar)];
end

if ismember(xVar,fieldnames(acq))
    xData=[acq.(xVar)];
end

%% Sort the data by the parameter given

[xData,inds]=sort(xData,'ascend');
atomdata=atomdata(inds);

%% Grab data

Nsum=zeros(length(atomdata),size(atomdata(1).RawImages,3));
for kk=1:length(atomdata)
    for jj=1:size(atomdata(1).RawImages,3)
       Nsum(kk,jj)=atomdata(kk).Raw.Total; 
    end
end

%% Linear Fit


if opts.FitLinear
    for jj=1:size(atomdata(1).RawImages,3)
       y=Nsum(:,jj);
       [pps{jj},pss{jj}]=polyfit(xData,y',1);
    end
end

%% Make Figure


% Create the name of the figure
[filepath,name,~]=fileparts(imgdir);

figDir=fullfile(imgdir,'figures');
if ~exist(figDir,'dir')
   mkdir(figDir); 
end

strs=strsplit(imgdir,filesep);
str=[strs{end-1} filesep strs{end}];

hF=figure('Name',['Raw Counts :'  str],...
    'units','pixels','color','w','Menubar','none','Resize','off');
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
ylabel('total counts');

hax.Position(4)=hax.Position(4)-20;

co=get(gca,'colororder');
coIN=brighten(co,.8);

Nbg=200*size(atomdata(1).RawImages,1)*size(atomdata(1).RawImages,2);
legStr={};
for nn=1:size(atomdata(1).RawImages,3)
   plot(xData,Nsum(:,nn)-0*Nbg,'o','color',co(nn,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(nn,:),'markeredgecolor',.5*co(nn,:),...
       'linewidth',2);
   ps(nn)=plot(xData,polyval(pps{nn},xData),'color',co(nn,:));
   legStr{nn}=[num2str(pps{nn}(1),'%.3e') ' counts/var'];
end

legend(ps,legStr,'location','northwest');


% Image directory folder string
t=uicontrol('style','text','string',str,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left');
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];

end

