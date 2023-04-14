function [hF2,outdata]=analyzeStripes2(ixondata,xVar,opts)

global ixon_imgdir
% maskname=fullfile('ixon_mask.mat');
% ixon_mask=load(maskname);
% ixon_mask=ixon_mask.BW;

doDebug=0;

if nargin==2
   opts=struct;
   opts.saveAnimation=1;
   opts.StartDelay=1;
   opts.MidDelay=0.2;
   opts.EndDelay=1;
   opts.ShimFit=0;
end
%% Sort the data by the parameter given
% Sort the data by the xVar and grab it.

params=[ixondata.Params];
xvals=[params.(xVar)];
[xvals,inds]=sort(xvals,'ascend');
ixondata=ixondata(inds);

%% Fitting Function
% Define the fitting function. It is a 2D gaussian who is modulated by a
% sine wave at an angle

gauss2dSine=@(A,xC,yC,sG,B,theta,L,phi,xx,yy) A*...
        exp(-((xx-xC).^2+(yy-yC).^2)/(2*sG^2)).*...
        (1+B*sin(2*pi/L*(cosd(theta)*xx+sind(theta)*yy)+phi));    

myfit=fittype(@(A,xC,yC,sG,B,theta,L,phi,xx,yy) ...
    gauss2dSine(A,xC,yC,sG,B,theta,L,phi,xx,yy),...
    'independent',{'xx','yy'},...
    'coefficients',{'A','xC','yC','sG','B','theta','L','phi'});

opt=fitoptions(myfit);

% Define fit options
% opt.Robust='bisquare';
opt.MaxIter=1000;
opt.MaxFunEvals=1000;
opt.TolFun=1E-8;


%% Make Filename
filename='ixon_animate_stripe2'; 
% Create the name of the figure
[filepath,name,~]=fileparts(ixon_imgdir);
figDir=fullfile(ixon_imgdir,'figures');
if ~exist(figDir,'dir')
   mkdir(figDir); 
end
% Make the figure name with the location
filename=fullfile(figDir,[filename '.gif']);

%% Construct initial guesses
% From the raw data, create the initial guesses for the stripe analysis.
%
%   X and Y center - gotten from center of mass
%   Stripe Angle - look for angle that generates largest contrast
%   Gaussian Raidus - look at second moment orthognal to the stripe angle
%   Wavelength - look at spacing of local maxima along summed counts anlong
%   the striep agnle
%   Phase - look at largerst correlator of data with a plane wave
%   Modulation depth - assume 50%

% Vectors for initial guess
xCG=zeros(length(ixondata),1);
yCG=zeros(length(ixondata),1);
thetaG=zeros(length(ixondata),1);
sG=zeros(length(ixondata),1);
BG=zeros(length(ixondata),1);
AG=zeros(length(ixondata),1);
phiG=zeros(length(ixondata),1);

fprintf('Constructing initial guesses ...');
t1=now;

% In debug mode, plot the guesses
if doDebug
    hf_guess=figure;
    clf
    hf_guess.Color='w';
    hf_guess.Position=[50 50 1000 400];
end

for kk=1:length(ixondata)    
    
    thetaVec=linspace(opts.Theta(1),opts.Theta(2),90);        
    
    % Get raw data
    Z2=ixondata(kk).Z;
    Z2(Z2<0)=0;
    x2=1:size(Z2,2);x2=x2';
    y2=1:size(Z2,1);y2=y2';    

    % Resize
    sc=0.1;
    Z2=imresize(Z2,sc);
    x2=imresize(x2,sc);
    y2=imresize(y2,sc);
    
    % Calculate cetner
    Zx=sum(Z2,1)'/sum(sum(Z2));    
    Zy=sum(Z2,2)/sum(sum(Z2));  
    Zx(Zx<0)=0;
    Zy(Zy<0)=0;
    
    xC=sum(Zx.*x2);
    yC=sum(Zy.*y2);     
    
    % Add guesses to external vector
    xCG(kk)=xC;
    yCG(kk)=yC;

    % Compute sum contrasts at many angles
    for jj=1:length(thetaVec)        
        Zrot=imrotate(Z2,thetaVec(jj));        
        Zsum=sum(Zrot,1);
        Zsum=smooth(Zsum,5);
        CC(jj)=sum(abs(diff(Zsum)).^2);    
    end
    
    % Find the angle which maximizes the contrast
    [~,ind]=max(CC);    
    theta=thetaVec(ind);
    thetaG(kk)=theta;

    % Find the cloud gaussian radius
    % Caclulate the second moment orthongal to the fringes
    Zrot=imrotate(Z2,theta,'crop');
    Zsum1=sum(Zrot,1)/sum(sum(Zrot));
    Zsum2=sum(Zrot,2)/sum(sum(Zrot));   
    Zsum2(Zsum2<0)=0;    
    yrC=sum(y2.*Zsum2);
    s=sqrt(sum(((y2-yrC).^2.*Zsum2)))*.75;
    
    sG(kk)=s;    
    
    % Find the guess wavelength
    % On the rotated data find the separation of local maxima
    ZsumSmooth=smooth(Zsum1,10);
    [yA,P]=islocalmax(ZsumSmooth,'MinSeparation',50*sc,...
        'MaxNumExtrema',4,'MinProminence',range(ZsumSmooth)*0.05);
    xA=diff(x2(yA));
    L=mean(xA);
    
    if isnan(L) || isinf(L) || isinf(-L) || L<50
       L=80; 
    end
    
    LG(kk)=L;
    
    % Find the initial phase of the data
    % Compute correlations of the image with a plane wave at the guess
    % angle
    phiVec=linspace(0,2*pi,50);
    [xx,yy]=meshgrid(x2,y2);
    Sphi=zeros(length(phiVec),1);
    for nn=1:length(phiVec)        
        plane_wave=sin(2*pi/L*(cosd(theta)*xx+sind(theta)*yy)+phiVec(nn));
        Sphi(nn)=sum(sum(plane_wave.*imgaussfilt(Z2,1)));     
    end
    [~,ind]=max(Sphi);
    phi=phiVec(ind);
    phiG(kk)=phi;

    % Guess the modulation strength
    B=max(P)/range(ZsumSmooth);
    B=0.5;
    BG(kk)=B;
    
    % Guess the gaussian envelope
    A=sum(sum(Z2))/(2*pi*s);
    A=max(max(Z2))/(2*(1+B));
    AG(kk)=A;
    
    % Construct the initial guess
    pG=[A,xC,yC,s,B,theta,L,phi];
    Zguess=gauss2dSine(pG(1),pG(2),pG(3),pG(4),pG(5),pG(6),pG(7),pG(8),xx,yy);
        
    % Plot the initial guess metrics    
    yL=[-1 1]*200.*sind(theta);
    xL=[-1 1]*200.*cosd(theta);
    
    yL2=[-1 1]*200.*sind(theta+90);
    xL2=[-1 1]*200.*cosd(theta+90);
    tt=linspace(0,2*pi,100);    
    %% Debug Plot for guess
    if doDebug
        % If in debug mode, plot the guess and some plots to show how it is
        % obtained.
        figure(hf_guess);
        clf
        colormap(purplemap);
        % Raw data
        subplot(231)
        imagesc(x2,y2,Z2);
        set(gca,'ydir','normal');
        axis equal tight
        hold on
        plot(xL+xC,yL+yC,'k-')    
        plot(xL2+xC,yL2+yC,'r-')    

        plot(s*cos(tt)+xC,s*sin(tt)+yC,'r-')
        colorbar
        
        % Show guess        
        subplot(232)   
        imagesc(x2,y2,Zguess);
        set(gca,'ydir','normal');
        axis equal tight
        colorbar

        % Guess residue
        subplot(233)  
        imagesc(x2,y2,Zguess-Z2);
        set(gca,'ydir','normal');
        axis equal tight
        colorbar
        drawnow;            
        
        % Contrast versus angle
        subplot(234)    
        plot(thetaVec,CC);
        xlabel('rotation angle (deg.)')
        ylabel('contrast (au)');
        xlim([min(thetaVec) max(thetaVec)]);
        hold on
        plot([1 1]*theta,get(gca,'YLim'),'r--');

        
        drawnow;    

        % Sum along non styriep direction local maxima
        subplot(235)    
        plot(x2,ZsumSmooth,'k-','linewidth',1);
        xlabel('rotated position (px)')
        ylabel('sum');
            hold on
        yP=yA.*ZsumSmooth;
        iP=[yP~=0];
        plot(x2(iP),yP(iP),'b*');
        xlim([min(x2) max(x2)]);
        plot(y2,Zsum2,'r-','linewidth',1);
        drawnow; 

        % Phase correlator
        subplot(236)    
        plot(phiVec/pi,Sphi);
        xlabel('initial phase (\pi)')
        ylabel('correlation');
        drawnow;  
%         keyboard
        waitforbuttonpress
    end
end
t2=now;
disp(['done (' num2str(round((t2-t1)*24*60*60),2) ' s)']);
%% Fit it

fRs={};
hF_live=figure;
hF_live.Color='w';
hF_live.Position=[100 100 1100 500];

co=get(gca,'colororder');

clf
colormap(purplemap);

ax1=subplot(4,4,[1 2 5 6 9 10 13 14]);    
hImg_raw=imagesc(ixondata(1).Z);
set(gca,'ydir','normal');
axis equal tight
% colorbar
    
% ax2=subplot(4,4,7);    
% hImg_fit=imagesc(ixondata(1).Z);
% set(gca,'ydir','normal');
% axis equal tight
% % colorbar

hold on
pFringe=plot(0,0,'-','color',co(1,:),'linewidth',2);
pPerp=plot(0,0,'-','color',co(5,:),'linewidth',2);     
pBar=plot(0,0,'-','color',co(1,:),'linewidth',2);     
pCirc=plot(0,0,'-','color',co(1,:),'linewidth',2);     


ax3=subplot(4,4,[11 15]);    
hImg_err=imagesc(ixondata(1).Z);
set(gca,'ydir','normal');
axis equal tight
% colorbar

ax4=subplot(4,4,[3 4 7 8]);    
pSum1_fit=plot(0,0,'r-','linewidth',2);
hold on
pSum1_data=plot(0,0,'-','color',co(1,:),'linewidth',2);
% set(gca,'ydir','normal');
xlabel('rotated position');
ylabel('sum counts');
%     
% ax5=subplot(4,4,[3 4]);    

pSum2_fit=plot(0,0,'r-','linewidth',1);
pSum2_data=plot(0,0,'-','color',co(5,:),'linewidth',1);
% set(gca,'ydir','normal');
% xlabel('position perpendicular');
% ylabel('sum counts');

ax6=subplot(4,4,[12 16]);
tbl=uitable('units','normalized','fontsize',8);
tbl.RowName={};
tbl.ColumnName={};
tbl.ColumnFormat={'char','char'};
tbl.ColumnWidth={120,80};

data={'Wavelength (px)', '';
    'Angle (deg.)','';
    'Phase (pi)', '';
    'Mod Depth','';
    'Xc', '';
    'Yc', '';
    'Amplitude', '';
    'Sigma (px)', ''};


tbl.Data=data;
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
sse=zeros(length(ixondata),1);


uicontrol('style','text','string','iXon, stripe','units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',12,'fontweight','bold',...
    'position',[2 2 100 20]);

vart=uicontrol('style','text','string','variable','units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',8,'fontweight','normal',...
    'position',[2 25 100 20]);

for kk=1:length(ixondata)
    fprintf(['Fitting ' num2str(kk) ' of ' num2str(length(ixondata)) ' ... ']);
    pG=[AG(kk),xCG(kk),yCG(kk),sG(kk),BG(kk),thetaG(kk),LG(kk),phiG(kk)];
    
    
    if kk>1 && isequal(xVar,'ExecutionDate')
        pG(end-2) = fout.theta;
        pG(end-1) = fout.L;
        pG(end) = fout.phi;    
    end
    
    opt.Lower = [0 0 0 0 0 pG(end-2)-2 pG(end-1)-3 pG(end)-1.5];
    opt.Upper = [AG(kk)*1.5 xCG(kk)+50 yCG(kk)+50 sG(kk)+60 BG(kk)+2000 pG(end-2)+2 pG(end-1)+3 pG(end)+1.5];


    
    opt.StartPoint=pG;

    
    % Grab the data
    Z=ixondata(kk).Z;
    Z=imgaussfilt(Z,1);
    x=1:size(Z,2);x=x';
    y=1:size(Z,1);
    
    % Create a meshgrid
    [xx,yy]=meshgrid(x,y);
    
    % Create the initial guess
    Zg=gauss2dSine(pG(1),pG(2),pG(3),pG(4),pG(5),pG(6),pG(7),pG(8),xx,yy); 
%     disp(pG)
    
    % Rescale the data for fitting
    sc=0.3;
    Zsc=imresize(Z,sc);
    xsc=imresize(x,sc);
    ysc=imresize(y,sc);
    [xxsc,yysc]=meshgrid(xsc,ysc);
    xxsc(Zsc<10)=[];
    yysc(Zsc<10)=[];
    Zsc(Zsc<10)=[];

    % Show the raw data and the initial guess
    set(hImg_raw,'XData',x,'YData',y,'CData',Z);
%     set(hImg_fit,'XData',x,'YData',y,'CData',Zg);
    set(hImg_err,'CData',Zg-Z);
    drawnow;   
    
    % Perform the fit
    [fout,gof,output]=fit([xxsc(:) yysc(:)],Zsc(:),myfit,opt);
    fRs{kk}=fout;
    sse(kk)=gof.sse;
    % Evalulate the fit
    Zf=feval(fout,xx,yy);
    disp(fout);
    % Overwrite the guess with the fit
    
%     set(hImg_fit,'CData',Zf);
    
    set(hImg_err,'CData',Zf-Z);
    set(ax1,'CLim',[0 fout.A*(1+fout.B)]);

%     set(ax2,'CLim',[0 fout.A]);
    set(ax3,'CLim',[-1 1]*fout.A*.5);

    set(pFringe,'XData',255+[0 1]*200*cosd(fout.theta),...
        'Ydata',255+[0 1]*200*sind(fout.theta));
    set(pPerp,'XData',255+[-1 1]*100*cosd(fout.theta+90),...
        'Ydata',255+[-1 1]*100*sind(fout.theta+90));
     set(pBar,'XData',255+[-50 50],...
        'Ydata',255*[1 1]);
    
    tt=linspace(0,fout.theta,100);
    set(pCirc,'XData',255+50*cosd(tt),...
        'YData',255+50*sind(tt));
    
    set(ax1,'XLim',[min(x) max(x)]);
%     set(ax2,'XLim',[min(x) max(x)]);
    set(ax3,'XLim',[min(x) max(x)]);

        
    % Show the sum counts along the stripe axis    
    set(pSum1_fit,'XData',x,'YData',sum(imrotate(Zf,fout.theta,'crop'),1));
    set(pSum1_data,'XData',x,'YData',sum(imrotate(Z,fout.theta,'crop'),1));
    set(ax4,'XLim',[min(1) max([max(x) max(y)])]);
    
    % Show the sum counts orthogonal to the stripe axis
    set(pSum2_fit,'XData',y,'YData',sum(imrotate(Zf,fout.theta,'crop'),2));
    set(pSum2_data,'XData',y,'YData',sum(imrotate(Z,fout.theta,'crop'),2));
%     set(ax5,'XLim',[min(y) max(y)]);
    
    % Summary table
    tbl.Data{1,2}=num2str(round(fout.L,3));
    tbl.Data{2,2}=num2str(round(fout.theta,3));    
    tbl.Data{3,2}=num2str(round(fout.phi/pi,3));
    tbl.Data{4,2}=num2str(round(fout.B,2));
    tbl.Data{5,2}=num2str(round(fout.xC,2));
    tbl.Data{6,2}=num2str(round(fout.yC,2));
    tbl.Data{7,2}=num2str(round(fout.A,2));
    tbl.Data{8,2}=num2str(round(fout.sG,2));
    drawnow;

    
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
    disp('done');
end

%% Process Fits
[cface1,cedge1] = ixoncolororder(1);

Ls=zeros(length(fRs),2);
thetas=zeros(length(fRs),2);
phis=zeros(length(fRs),2);

for kk=1:length(fRs)
    ci=confint(fRs{kk});
    
    % Get fit value
    Ls(kk,1)=fRs{kk}.L; 
    thetas(kk,1)=fRs{kk}.theta;
    phis(kk,1)=fRs{kk}.phi;

    % Get the 95% confidence interval
    Ls(kk,2)=(ci(2,7)-ci(1,7))/2;   
    thetas(kk,2)=(ci(2,6)-ci(1,6))/2;    
    phis(kk,2)=(ci(2,8)-ci(1,8))/2;       

end

% Summarize the results
hF2=figure;
hF2.Position=[100 100 900 250];
hF2.Color='w';
clf;

uicontrol('style','text','string','iXon, stripe','units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',12,'fontweight','bold',...
    'position',[2 2 100 20]);

% Folder directory
strs=strsplit(ixon_imgdir,filesep);
str=[strs{end-1} filesep strs{end}];
% Image directory folder string
t=uicontrol('style','text','string',str,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',6);
t.Position(4)=t.Extent(4);
t.Position(3)=hF2.Position(3);
t.Position(1:2)=[5 hF2.Position(4)-t.Position(4)];


% Wavelength
subplot(131);


errorbar(xvals,Ls(:,1),Ls(:,2),'marker','o',...
    'MarkerFacecolor',cface1,'markeredgecolor',cedge1,'linestyle','none',...
    'linewidth',1.5,'color',cedge1);
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
ylabel('wavelength (px)');
grid on

if isequal(xVar,'ExecutionDate')
    datetick('x');
    xlabel('time','fontsize',10);
end
xlim([min(xvals) max(xvals)]);

% Angle
axb2=subplot(132);
errorbar(xvals,thetas(:,1),thetas(:,2),'marker','o',...
    'MarkerFacecolor',cface1,'markeredgecolor',cedge1,'linestyle','none',...
    'linewidth',1.5,'color',cedge1);
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
ylabel('angle (deg.)');
grid on
if isequal(xVar,'ExecutionDate')
    datetick('x');
    xlabel('time','fontsize',10);
end
xlim([min(xvals) max(xvals)]);

% Phase
subplot(133);
errorbar(xvals,phis(:,1)/pi,phis(:,2)/pi,'marker','o',...
    'MarkerFacecolor',cface1,'markeredgecolor',cedge1,'linestyle','none',...
    'linewidth',1.5,'color',cedge1);
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
ylabel('phase (\pi)');
grid on
if isequal(xVar,'ExecutionDate')
    datetick('x');
    xlabel('time','fontsize',10);
end
xlim([min(xvals) max(xvals)]);

% Error
% subplot(144);
% plot(xvals,sse,'marker','o',...
%     'MarkerFacecolor',cface1,'markeredgecolor',cedge1,'linestyle','none',...
%     'linewidth',1.5,'color',cedge1);
% xlabel(xVar);
% ylabel('sse (au)');

if opts.ShimFit
    % This fitting algotrithm is specifically for when a shim varied and
    % the other one is held constant
%     wave_func=@(A,x0,x1,x) A*sqrt(1/((x-x0).^2+x1.^2));
%     theta_func=@(A,x0,x1,x) atan2(x-x0,x1);
%     
    switch xVar
        case 'xshimd'
            theta_c=59.8;            
        case 'yshimd'
            theta_c=-31.3;            
        otherwise
            warning('uh oh, go ask cora what went wrong');
            theta_c=-45;        
    end
    
    % Find angles close by to the critical one
    inds=[abs(thetas(:,1)-theta_c)<10];
    
    [~,i0]=min(xvals);
    b0=thetas(i0,1);
   
    
    xp=xvals(inds)';
    yp=thetas(inds,1);
    
    
    
    m0=range(yp)./range(xp);
    
    myfit2=fittype('m*x+b','independent','x','coefficients',{'m','b'});
    opt2=fitoptions(myfit2);
    opt2.StartPoint=[m0 b0];
    
    fout_theta=fit(xp,yp,myfit2,opt2);
    
    axes(axb2)
    hold on
    tt=linspace(min(xp),max(xp),10);
    pFit=plot(tt,feval(fout_theta,tt),'k-','linewidth',1);
    
    pstr=['(' num2str(round(fout_theta.m,3)) ')x + ' num2str(round(fout_theta.b,3))];
    
    legend(pFit,pstr,'location','best');
    
end

%% Prepare output data
outdata=struct;
outdata.xVar=xVar;
outdata.xvals=xvals;
outdata.Wavelength=Ls;
outdata.Theta=thetas;
outdata.Phi=phis;
end
