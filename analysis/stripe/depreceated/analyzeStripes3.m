function stripes = analyzeStripes3(ixondata,xVar,opts)

global ixon_imgdir

if nargin==2
   opts=struct;
   opts.saveAnimation=1;
   opts.StartDelay=1;
   opts.MidDelay=0.05;
   opts.EndDelay=1;
end

ind = 1;

%% Sort the data by the parameter given
% Sort the data by the xVar and grab it.

params=[ixondata.Params];
xvals=[params.(xVar)];
[xvals,inds]=sort(xvals,'ascend');
ixondata=ixondata(inds);

%% Make Filename
filename='ixon_animate_stripe'; 
% Create the name of the figure
[filepath,name,~]=fileparts(ixon_imgdir);
figDir=fullfile(ixon_imgdir,'figures');
if ~exist(figDir,'dir')
   mkdir(figDir); 
end
% Make the figure name with the location
filename=fullfile(figDir,[filename '.gif']);

%% Initialize Figure

% Make the figure
hF_live=figure;
hF_live.Color='w';
hF_live.Position=[100 500 1100 500];
co=get(gca,'colororder');
clf
colormap(purplemap);

% Image Plot
ax1=subplot(4,4,[1 2 5 6 9 10 13 14]);    
hImg_raw=imagesc(ixondata(1).X,ixondata(1).Y,ixondata(1).Z(:,:,ind));
set(gca,'ydir','normal');
axis equal tight  
hold on
pFringe=plot(0,0,'-','color',co(1,:),'linewidth',2);
pPerp=plot(0,0,'-','color',co(5,:),'linewidth',2);     
pBar=plot(0,0,'-','color',co(1,:),'linewidth',2);     
pCirc=plot(0,0,'-','color',co(1,:),'linewidth',2);     

% Residue plot
ax3=subplot(4,4,[11 15]);    
hImg_err=imagesc(ixondata(1).X,ixondata(1).Y,ixondata(1).Z(:,:,ind));
set(gca,'ydir','normal');
axis equal tight
% colorbar

% Sum Plots
ax4=subplot(4,4,[3 4 7 8]);    
pSum2_fit=plot(0,0,'k--','linewidth',1);
hold on
pSum2_data=plot(0,0,'-','color',co(5,:),'linewidth',1);

pSum1_fit=plot(0,0,'-','linewidth',2,'color',co(2,:));
hold on
pSum1_data=plot(0,0,'-','color',co(1,:),'linewidth',1);
xlabel('rotated position (px)');
ylabel('sum counts');

str = {'fringe','fringe fit','perp','perp fit'};
legend([pSum1_data, pSum1_fit, pSum2_data, pSum2_fit],str,'fontsize',8,...
    'location','northeast')
text(.01,.98,'projected sum counts','units','normalized',...
    'verticalalignment','top','fontsize',8);
% Output Table
ax6=subplot(4,4,[12 16]);
tbl=uitable('units','normalized','fontsize',8);
tbl.RowName={};
tbl.ColumnName={};
tbl.ColumnFormat={'char','char'};
tbl.ColumnWidth={120,80};
tbl_labels = {'Amplitude', '';
    'xC (px)', '';
    'yC (px)', '';
    'Sigma1 (px)', '';
    'Sigma2 (px)', '';
    'Mod Depth', '';
    'Rotation (deg)', '';
    'Wavelength (px)', '';
    'Phase (2*pi)', '';};
tbl.Data=tbl_labels;
tbl.Position(3:4)=tbl.Extent(3:4);
tbl.Position(1:2)=ax6.Position(1:2);
delete(ax6);    
      

% Folder directory
strs=strsplit(ixon_imgdir,filesep);
str=[strs{end-1} filesep strs{end}];
% Image directory folder string
t=uicontrol('style','text','string',str,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',6);
t.Position(4)=t.Extent(4);
t.Position(3)=hF_live.Position(3);
t.Position(1:2)=[5 hF_live.Position(4)-t.Position(4)];

uicontrol('style','text','string','iXon, stripe','units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',12,'fontweight','bold',...
    'position',[2 2 100 20]);

vart=uicontrol('style','text','string','variable','units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',8,'fontweight','normal',...
    'position',[125 2 400 20]);
%% Iterate fit and update graphics

clear stripes
for kk=1:length(ixondata)
    fprintf([num2str(kk) '/' num2str(length(ixondata)) ' ']);
   
    % Grab the data
    z = ixondata(kk).Z(:,:,ind);
    x = ixondata(kk).X;
    y = ixondata(kk).Y;
    [xx,yy] = meshgrid(x,y);


    % Show the raw data and the initial guess       

    % Fit the stripes
    stripes(kk)=ixon_fitStripe(x',y',z,opts);
    theta = stripes(kk).theta;

    % Evalulate the fit
    Zfit=feval(stripes(kk).Fit,xx,yy);

    % Overwrite the guess with the fit   
    set(hImg_raw,'XData',x,'YData',y,'CData',z);
    set(hImg_err,'CData',Zfit-z);
    set(ax1,'CLim',opts.CLim);

    set(pFringe,'XData',255+[0 1]*200*cosd(theta),...
        'Ydata',255+[0 1]*200*sind(theta));
    set(pPerp,'XData',255+[-1 1]*100*cosd(theta+90),...
        'Ydata',255+[-1 1]*100*sind(theta+90));
     set(pBar,'XData',255+[-50 50],...
        'Ydata',255*[1 1]);
    
    tt=linspace(0,theta,100);
    set(pCirc,'XData',255+50*cosd(tt),...
        'YData',255+50*sind(tt));
    
    set(ax1,'XLim',[min(x) max(x)]);
    set(ax3,'XLim',[min(x) max(x)]);
        
    % Show the sum counts along the stripe axis    
    set(pSum1_fit,'XData',x,'YData',sum(imrotate(Zfit,theta,'crop'),1));
    set(pSum1_data,'XData',x,'YData',sum(imrotate(z,theta,'crop'),1));
    set(ax4,'XLim',[min(1) max([max(x) max(y)])]);
    
    % Show the sum counts orthogonal to the stripe axis
    set(pSum2_fit,'XData',y,'YData',sum(imrotate(Zfit,theta,'crop'),2));
    set(pSum2_data,'XData',y,'YData',sum(imrotate(z,theta,'crop'),2));

    set(ax4,'YLim',[0 ax4.YLim(2)]);

    % Summary table
    tbl.Data{1,2}=num2str(round(stripes(kk).A,2));
    tbl.Data{2,2}=num2str(round(stripes(kk).xC,2));
    tbl.Data{3,2}=num2str(round(stripes(kk).yC,2));
    tbl.Data{4,2}=num2str(round(stripes(kk).s1,2));
    tbl.Data{5,2}=num2str(round(stripes(kk).s2,2));
    tbl.Data{6,2}=num2str(round(stripes(kk).B,2));
    tbl.Data{7,2}=num2str(round(stripes(kk).theta,3));       
    tbl.Data{8,2}=num2str(round(stripes(kk).L,3));
    tbl.Data{9,2}=num2str(round(stripes(kk).phi/(2*pi),3));
    drawnow;    
    
    if isequal(xVar,'ExecutionDate')
        vart.String=[xVar ': ' datestr(xvals(kk),'YYYY-mm-DD_HH-MM-SS')];          % Variable string
    else
        vart.String=[xVar ': ' num2str(xvals(kk))];          % Variable string
    end
    
    if opts.saveAnimation    
         % Write the image data
        frame = getframe(hF_live);
        im = frame2im(frame);
        [A,map] = rgb2ind(im,256);           
        switch kk
            case 1
                imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',opts.StartDelay);
            case length(ixondata)
                imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',opts.EndDelay);
            otherwise
                imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',opts.MidDelay);
        end
    end
end

%% Unwrap the phase


% %% Prepare output data
% outdata=struct;
% outdata.xVar=xVar;
% outdata.xvals=xvals;
% outdata.Wavelength=Ls;
% outdata.Theta=thetas;
% outdata.Phi=phis;
end
