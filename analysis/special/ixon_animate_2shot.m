function ixon_animate_2shot(ixondata,xVar,opts)
clim=opts.CLim;
global ixon_imgdir


%% Sort Data
% Sort the data by xVar in either ascending or descending order
direction=opts.Order;

params=[ixondata.Params];
xvals=[params.(xVar)];
if isequal(direction,'ascend')
    [xvals,inds]=sort(xvals,'ascend');
else
    [xvals,inds]=sort(xvals,'descend');
end
ixondata=ixondata(inds);
%% Auto Clim

if isequal(clim,'auto')
  cL0=0;
  cH0=1;
    for kk=1:size(ixondata,3)
       zH = max(ixondata(kk).(opts.Source),[],'all');
       cH0 = max([zH cH0]);    
    end
    cH0 = cH0*1.2;
    clim=[cL0 cH0];
    
end

%% Get First Data Point

% grab initial data
Z=ixondata(1).(opts.Source)(:,:,1);
Y=ixondata(1).Y;
X=ixondata(1).X;

%% Animation Settings
startDelay=opts.StartDelay;   % First picture hold time
midDelay=opts.MidDelay;   % Middle pictures hold time
endDelay=opts.EndDelay;     % End picture hold time

strs=strsplit(ixon_imgdir,filesep);
str=[strs{end-1} filesep strs{end}];

%% Make Filename
if ~isfield(opts,'filename')
    filename='ixon_animate_2shot'; 
else
    filename= opts.filename;
end

% Create the name of the figure
[filepath,name,~]=fileparts(ixon_imgdir);

figDir=fullfile(ixon_imgdir,'figures');
if ~exist(figDir,'dir')
   mkdir(figDir); 
end

% Make the figure name with the location
filename=fullfile(figDir,[filename '.gif']);

%% Make Figure

hF=figure('Name',[str ' : Ixon Animate'],...
    'units','pixels','color','w','Menubar','none','Resize','off',...
    'WindowStyle','modal');
hF.Position(1)=10;
hF.Position(2)=5;
hF.Position(3)=1200;
hF.Position(4)=550;
drawnow;

% Image directory folder string
t=uicontrol('style','text','string',str,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left');
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];

% First axis
ax1 = subplot(121);
set(ax1,'box','on','xgrid','on','ygrid','on','ydir','normal');
hImg1=imagesc(X,Y,Z);
axis equal tight
caxis(clim);
cend = [0.6 0 .5];
colormap([linspace(1,cend(1),1000)' linspace(1,cend(2),1000)' linspace(1,cend(3),1000)'])
hold on
colorbar
set(gca,'units','pixels','box','on','linewidth',1,'ydir','normal');
% title('image 1');
% Text label for variable name
t1_a=text(.01,.01,'hi','units','normalized','fontsize',12,'color','r',...
    'interpreter','none','verticalalignment','bottom',...
    'fontweight','bold');
t1_b=text(.01,.99,'hi','units','normalized','fontsize',12,'color','r',...
    'interpreter','none','verticalalignment','top',...
    'fontweight','bold');


% Second axis
ax2 = subplot(122);
set(ax2,'box','on','xgrid','on','ygrid','on','ydir','normal');
hImg2=imagesc(X,Y,Z);
axis equal tight
caxis(clim);
cend = [0.6 0 .5];
colormap([linspace(1,cend(1),1000)' linspace(1,cend(2),1000)' linspace(1,cend(3),1000)']);
hold on
colorbar
set(gca,'units','pixels','box','on','linewidth',1,'ydir','normal');
% title('image 2');
t2_a=text(.01,.01,'hi','units','normalized','fontsize',12,'color','r',...
    'interpreter','none','verticalalignment','bottom',...
    'fontweight','bold');
t2_b=text(.01,.99,'hi','units','normalized','fontsize',12,'color','r',...
    'interpreter','none','verticalalignment','top',...
    'fontweight','bold');


%% Animate
disp('ANIMATING');
for kk=1:length(ixondata)   % Iterate over all unique xvalues
    
    set(t1_b,'String',[ixondata(kk).Name ' (1)']);
    set(t2_b,'String',[ixondata(kk).Name  ' (2)']);

    %%%% Update the graphics
    
    if isequal(xVar,'ExecutionDate')
        t1_a.String=[xVar ': ' datestr(xvals(kk),'YYYY-mm-DD_HH-MM-SS') ' (1)'];          % Variable string
        t2_a.String=[xVar ': ' datestr(xvals(kk),'YYYY-mm-DD_HH-MM-SS') ' (2)'];          % Variable string

    else
        t1_a.String=[xVar ': ' num2str(xvals(kk))];          % Variable string
        t2_a.String=[xVar ': ' num2str(xvals(kk))];          % Variable string

    end
    
    set(hImg1,'XData',X,'YData',Y,'CData',ixondata(kk).(opts.Source)(:,:,1));  % Image data
    set(hImg2,'XData',X,'YData',Y,'CData',ixondata(kk).(opts.Source)(:,:,2));  % Image data

    drawnow % update graphcis
    
    
    % Write the image data
    frame = getframe(hF);
    im = frame2im(frame);
    [A,map] = rgb2ind(im,256);           

    if kk == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',startDelay);
    else
        if kk==length(ixondata)
            imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',endDelay);
        else
            imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',midDelay);
        end
    end

end
close;
    disp('done');
end

