function ixon_animate(ixondata,xVar,opts)
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


%% Average with unique x variable

% Find and sor the unique values
uxvals=unique(xvals);

if isequal(direction,'ascend')
    uxvals=sort(uxvals,'ascend');
else
    uxvals=sort(uxvals,'descend');
end

Zall=zeros(size(ixondata(1).Z,1),size(ixondata(1).Z,2),length(uxvals));

for kk=1:length(uxvals) % Iterate over unique x values    
    % Find the indeces which have this unique value
    inds=find(uxvals(kk)==xvals);    
    for ii=1:length(inds)
        ind=inds(ii);
        Z=ixondata(ind).Z;
        Zall(:,:,kk)=Zall(:,:,kk)+Z;        
    end        
    Zall(:,:,kk)=Zall(:,:,kk)/length(inds);   
end

%% Auto Clim

if isequal(clim,'auto')
  cL0=0;
  cH0=1;
    for kk=1:size(Zall,3)
        Z=Zall(:,:,kk);
        N0=size(Z,1)*size(Z,2);
        dN=round(N0*.01); % pad by 1% floor
        call=sort(Z(:));
        cL=call(dN);
        cH=call(end);     
        cH0=max([cH0 cH]);
        cL0=min([cL0 cL]);
    end
    clim=[cL0 cH0];
end

%% Animation Settings
startDelay=opts.StartDelay;   % First picture hold time
midDelay=opts.MidDelay;   % Middle pictures hold time
endDelay=opts.EndDelay;     % End picture hold time

strs=strsplit(ixon_imgdir,filesep);
str=[strs{end-1} filesep strs{end}];

%% Make Filename
filename='ixon_animate'; 

% Create the name of the figure
[filepath,name,~]=fileparts(ixon_imgdir);

figDir=fullfile(ixon_imgdir,'figures');
if ~exist(figDir,'dir')
   mkdir(figDir); 
end

% Make the figure name with the location
filename=fullfile(figDir,[filename '.gif']);

%% Make Figure

% grab initial data
Z=ixondata(1).Z;
Y=1:size(Z,1);
X=1:size(Z,2);


hF=figure('Name',[str ' : Ixon Animate'],...
    'units','pixels','color','w','Menubar','none','Resize','off',...
    'WindowStyle','modal');
hF.Position(1)=10;
hF.Position(2)=5;
hF.Position(3)=600;
hF.Position(4)=550;
drawnow;

% Image directory folder string
t=uicontrol('style','text','string',str,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left');
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];

% Axes for data
hAxImg=axes('parent',hF,'units','pixels','Box','on','XGrid','on',...
    'YGrid','on','YDir','normal','XAxisLocation','bottom');
hAxImg.Position(2)=50;
hAxImg.Position(1)=60;
hAxImg.Position(3)=hAxImg.Position(4);
hAxImg.Position(3:4)=hF.Position(3:4)-[80 80];
drawnow;


% Text label for variable name
t=text(2,2,'hi','units','pixels','fontsize',12,'color','r',...
    'interpreter','none','verticalalignment','bottom',...
    'fontweight','bold');

%colormap(purplemap)
cend = [0.6 0 .5];
colormap([linspace(1,cend(1),1000)' linspace(1,cend(2),1000)' linspace(1,cend(3),1000)'])

hold on
hImg=imagesc(X,Y,Z);
axis equal tight
caxis(clim);
co=get(gca,'colororder');

hold on
colorbar

set(gca,'units','pixels','box','on','linewidth',2,'ydir','normal');

% Add ROI
for kk=1:size(ixondata(1).ROI,1)
    ROI=ixondata(1).ROI(kk,:);
    x0=ROI(1);
    y0=ROI(3);
    H=ROI(4)-ROI(3);
    W=ROI(2)-ROI(1);
    pROI=rectangle('position',[x0 y0 W H],'edgecolor',co(kk,:),'linewidth',2);
end

drawnow;

% keyboard

%% Animate
disp('ANIMATING');
for kk=1:length(uxvals)   % Iterate over all unique xvalues
    
    %%%% Update the graphics
    
    if isequal(xVar,'ExecutionDate')
            t.String=[xVar ': ' datestr(uxvals(kk),'YYYY-mm-DD_HH-MM-SS')];          % Variable string
    else
            t.String=[xVar ': ' num2str(uxvals(kk))];          % Variable string
    end


    set(hImg,'XData',X,'YData',Y,'CData',Zall(:,:,kk));  % Image data
    set(gca,'XDir','normal','YDir','normal');
    
    drawnow % update graphcis
    
    
    % Write the image data
    frame = getframe(hF);
    im = frame2im(frame);
    [A,map] = rgb2ind(im,256);           

    if kk == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',startDelay);
    else
        if kk==length(uxvals)
            imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',endDelay);
        else
            imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',midDelay);
        end
    end

end
close;
    disp('done');
end

