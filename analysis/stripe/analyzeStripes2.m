function [hF2,outdata]=analyzeStripes2(ixondata,xVar,opts)
tic
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

% [xvals,inds]=sort(xvals,'descend');

ixondata=ixondata(inds);

%% Fitting Function
% Define the fitting function. It is a 2D gaussian who is modulated by a
% sine wave at an angle

phase_map = @(L,theta,phi,xx,yy) ...
    2*pi/L*(cosd(theta)*xx+sind(theta)*yy)+phi;

gauss2dSine=@(A,xC,yC,sG,B,theta,L,phi,xx,yy) A*...
        exp(-((xx-xC).^2+(yy-yC).^2)/(2*sG^2)).*...
        (1+B*sin(phase_map(L,theta,phi,xx,yy)));    

myfit=fittype(@(A,xC,yC,sG,B,theta,L,phi,xx,yy) ...
    gauss2dSine(A,xC,yC,sG,B,theta,L,phi,xx,yy),...
    'independent',{'xx','yy'},...
    'coefficients',{'A','xC','yC','sG','B','theta','L','phi'});
opt=fitoptions(myfit);

% mod_func = @(B,theta,L,phi,duty,xx,yy) ...
%     (1-B*(0.5*square(phase_map(L,theta,phi,xx,yy)+pi/2,duty)+.5));
% gauss2dsquare=@(A,xC,yC,sG,B,theta,L,phi,duty,xx,yy) A*...
%         exp(-((xx-xC).^2+(yy-yC).^2)/(2*sG^2)).*...
%         mod_func(B,theta,L,phi,duty,xx,yy);    

% 
% myfit2=fittype(@(A,xC,yC,sG,B,theta,L,phi,duty,xx,yy) ...
%     gauss2dsquare(A,xC,yC,sG,B,theta,L,phi,duty,xx,yy),...
%     'independent',{'xx','yy'},...
%     'coefficients',{'A','xC','yC','sG','B','theta','L','phi','duty'});
% opt2=fitoptions(myfit);

% https://en.wikipedia.org/wiki/Gaussian_function
% But we add a minus sign to make it counter clockwise angle
% When theta=0 s1 is on the x axis
a = @(theta,s1,s2) cosd(theta)^2/(2*s1^2) + sind(theta)^2/(2*s2^2);
b = @(theta,s1,s2) -sind(2*theta)/(4*s1^2) + sind(2*theta)/(4*s2^2);
c = @(theta,s1,s2) sind(theta)^2/(2*s1^2) + cosd(theta)^2/(2*s2^2);

gauss2dSineRot = @(A,xC,yC,s1,s2,B,theta,L,phi,xx,yy) ...
    A.*exp(-(...
    a(theta,s1,s2)*(xx-xC).^2 + ...
    2*b(theta,s1,s2)*(xx-xC).*(yy-yC)+ ...
    c(theta,s1,s2)*(yy-yC).^2)).*...
    (1+B*sin(2*pi/L*(cosd(theta)*xx+sind(theta)*yy)+phi));

myfit2=fittype(@(A,xC,yC,s1,s2,B,theta,L,phi,xx,yy) ...
    gauss2dSineRot(A,xC,yC,s1,s2,B,theta,L,phi,xx,yy),...
    'independent',{'xx','yy'},...
    'coefficients',{'A','xC','yC','s1','s2','B','theta','L','phi'});
opt2=fitoptions(myfit);

             

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

%% Fit it

fRs={};
hF_live=figure;
hF_live.Color='w';
hF_live.Position=[100 500 1100 500];

co=get(gca,'colororder');

clf
colormap(purplemap);

ax1=subplot(4,4,[1 2 5 6 9 10 13 14]);    
hImg_raw=imagesc(ixondata(1).X,ixondata(1).Y,ixondata(1).Z);
set(gca,'ydir','normal');
axis equal tight
   

hold on
pFringe=plot(0,0,'-','color',co(1,:),'linewidth',2);
pPerp=plot(0,0,'-','color',co(5,:),'linewidth',2);     
pBar=plot(0,0,'-','color',co(1,:),'linewidth',2);     
pCirc=plot(0,0,'-','color',co(1,:),'linewidth',2);     


ax3=subplot(4,4,[11 15]);    
hImg_err=imagesc(ixondata(1).X,ixondata(1).Y,ixondata(1).Z);
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

data={'Wavelength (px)', '';
    'Angle (deg.)','';
    'Phase (2*pi)', '';
    'Mod Depth','';
    'Xc', '';
    'Yc', '';
    'Amplitude', '';
    'Sigma1 (px)', '';
    'Sigma2 (px)', ''};

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
    'position',[125 2 400 20]);
toc
%%

% Initialize outputs (value,error)
% A,xC,yC,s1,s2,B,theta,L,phi
As = zeros(length(ixondata),2);
xCs = zeros(length(ixondata),2);
yCs = zeros(length(ixondata),2);
s1s = zeros(length(ixondata),2);
s2s = zeros(length(ixondata),2);
thetas=zeros(length(ixondata),2);
Ls=zeros(length(ixondata),2);
phis=zeros(length(ixondata),2);
R2s = zeros(length(ixondata),1);
sse = zeros(length(ixondata),1);


for kk=1:length(ixondata)
    fprintf([num2str(kk) '/' num2str(length(ixondata)) ' ']);
      
    % Grab the data
    Z=ixondata(kk).Z;
    Z=imgaussfilt(Z,1);
    x=ixondata(kk).X;x=x';y =ixondata(kk).Y;y=y'; 
    
    % Create a 2D meshgrid
    [xx,yy]=meshgrid(x,y);

    % Make the initial guess
    tic;G=calculateInitialGuess(x,y,Z,opts);t=toc;
    opt2.StartPoint = G;
    % fprintf(['(' num2str(t,2) 's) ']);

    % Assign reasonable upper and lower bounds
    % A,xC,yC,s1,s2,B,theta,L,phi
    opt2.Upper = [G(1)*2 512 512 250 250 1  360 512 inf];
    opt2.Lower = [0        0   0   0   0 0 -360   0 -inf];

    % If secondary data point, choose an initial phase which is within
    % +-pi of previous.
    if kk>1         
        phi0 = opt2.StartPoint(end);        
        phivec = phi0+(-100:1:100)*2*pi;
        phiold = fout.phi;         
        [val,ind]=min(abs(phivec-phiold));        
        phinew = phivec(ind);
        opt2.StartPoint(end) = phinew;  
    end                     
     
    Zg=gauss2dSineRot(G(1),G(2),G(3),G(4),G(5),G(6),G(7),G(8),G(9),xx,yy);

    % Rescale the data for fitting
    sc=0.3;
    Zsc=imresize(Z,sc);xsc=imresize(x,sc);ysc=imresize(y,sc);
    [xxsc,yysc]=meshgrid(xsc,ysc);
    
    if isfield(opts,'Threshhold') && ~isnan(opts.Threshhold)
        lvl = opts.Threshhold;
        xxsc(Zsc<lvl)=[];yysc(Zsc<lvl)=[];Zsc(Zsc<lvl)=[];  
    else
        xxsc(Zsc<=0)=[];yysc(Zsc<=0)=[];Zsc(Zsc<=0)=[];        
    end
    
    
    if isfield(opts,'ROI') && ~sum(isnan(opts.ROI))
        ROI = opts.ROI;        
        ii = xxsc<ROI(1);
        xxsc(ii)=[];yysc(ii)=[];Zsc(ii)=[];          
        i2 = xxsc>ROI(2);
        xxsc(ii)=[];yysc(ii)=[];Zsc(ii)=[];  
        i3 = yysc<ROI(3);
        xxsc(ii)=[];yysc(ii)=[];Zsc(ii)=[];  
        i4 = yysc>ROI(4);
        xxsc(ii)=[];yysc(ii)=[];Zsc(ii)=[];  
    end

    % Show the raw data and the initial guess
    set(hImg_raw,'XData',x,'YData',y,'CData',Z);
%     set(hImg_fit,'XData',x,'YData',y,'CData',Zg);
    set(hImg_err,'CData',Zg-Z);
    drawnow;   
    % Perform the fit
    % fprintf('fitting ...');
    fprintf('(A,xC,yC,s1,s2,B,theta,L,phi):');
    s = ['(' ...
        num2str(round(opt2.StartPoint(1))) ',' ...
        num2str(round(opt2.StartPoint(2))) ',' ...
        num2str(round(opt2.StartPoint(3))) ',' ...
        num2str(round(opt2.StartPoint(4))) ',' ...
        num2str(round(opt2.StartPoint(5))) ',' ...
        num2str(round(opt2.StartPoint(6))) ',' ...
        num2str(round(opt2.StartPoint(7))) ',' ...
        num2str(round(opt2.StartPoint(8))) ',' ...
        num2str(round(opt2.StartPoint(9),2)) ')'];
    fprintf(s);
    fprintf('->');
    tic;
    [fout,gof,output]=fit([xxsc(:) yysc(:)],Zsc(:),myfit2,opt2);    
    t=toc;    

    s = ['(' ...
        num2str(round(fout.A)) ',' ...
        num2str(round(fout.xC)) ',' ...
        num2str(round(fout.xC)) ',' ...
        num2str(round(fout.s1)) ',' ...
        num2str(round(fout.s2)) ',' ...
        num2str(round(fout.B)) ',' ...
        num2str(round(fout.theta)) ',' ...
        num2str(round(fout.L)) ',' ...
        num2str(round(fout.phi,2)) ') '];
    fprintf(s)
    disp(['(' num2str(t,2) 's) ']);                
     
    

    % Make sure second data point is within +-pi
   if kk==2         
       if abs(fout.phi-phiold)>pi
           warning('Fitting again to fix the phase')
            G = [fout.A fout.xC fout.yC fout.s1 fout.s2 fout.B fout.theta fout.L fout.phi];
            if fout.phi>phiold 
                G(end) = G(end)-2*pi;
            else 
                G(end) = G(end)+2*pi;
            end
            opt2.StartPoint = G;
            tic;
            [fout,gof,output]=fit([xxsc(:) yysc(:)],Zsc(:),myfit2,opt2);    
            t=toc;    
       end     
   end   

   % Make sure that large phase jumps conserve the slope sign dphi/dvar
    if kk>2         
        try
       if abs(fout.phi-phiold)>.8*pi     
           warning('large phase jump detected, attempting to preserver the slope')
            thisX = xvals(kk);      % Current X value

            % Get the unique xvalues
            [ux,inds]=unique(xvals);
            uphi = phis(inds);
            iMe = find(ux==thisX,1);

            % Find the previous two phases and compute slope
            phi1 = uphi(iMe-1);x1 = ux(iMe-1);
            phi2 = uphi(iMe-2);x2 = ux(iMe-2);         
            dPhi = (phi1-phi2)/(x1-x2);

            dXme = sign(thisX-x1);

            % Create a list of new possible phis to fit
            phiVec = fout.phi + [-100:1:100]*2*pi;

            % Find the first phi is that is larger than previous
            ind=find(phiVec>phis(kk-1),1);

            phiLarger = phiVec(ind);
            phiSmaller = phiVec(ind-1);

            G = [fout.A fout.xC fout.yC fout.s1 fout.s2 fout.B fout.theta fout.L fout.phi];

            if sign(dPhi*dXme)==1
                G(end) = phiLarger;
            else
                G(end) = phiSmaller;
            end
       
            opt2.StartPoint = G;
            tic;
            [fout,gof,output]=fit([xxsc(:) yysc(:)],Zsc(:),myfit2,opt2);    
            t=toc;    
       end   
        catch ME
            warning(ME);
        end
   end  


    % Process output
    ci=confint(fout);
    fRs{kk}=fout;
    As(kk,:) = [fout.A (ci(2,1)-ci(1,1))*.5];
    xCs(kk,:) = [fout.xC (ci(2,2)-ci(1,2))*.5];
    yCs(kk,:) =  [fout.yC (ci(2,3)-ci(1,3))*.5];
    s1s(kk,:) =  [fout.s1 (ci(2,4)-ci(1,4))*.5];
    s2s(kk,:) = [fout.s2 (ci(2,5)-ci(1,5))*.5];
    thetas(kk,:)=[fout.theta (ci(2,6)-ci(1,6))*.5];
    Ls(kk,:)=[fout.L (ci(2,7)-ci(1,7))*.5];
    phis(kk,:)=[fout.phi (ci(2,8)-ci(1,8))*.5];
    
    R2s(kk) = gof.rsquare;
    sse(kk) = gof.sse;


    % Evalulate the fit
    Zf=feval(fout,xx,yy);

    % Overwrite the guess with the fit   
    set(hImg_err,'CData',Zf-Z);
    set(ax1,'CLim',opts.CLim);


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
    tbl.Data{3,2}=num2str(round(fout.phi/(2*pi),3));
    tbl.Data{4,2}=num2str(round(fout.B,2));
    tbl.Data{5,2}=num2str(round(fout.xC,2));
    tbl.Data{6,2}=num2str(round(fout.yC,2));
    tbl.Data{7,2}=num2str(round(fout.A,2));
%     tbl.Data{8,2}=num2str(round(fout.sG,2));

    tbl.Data{8,2}=num2str(round(fout.s1,2));
    tbl.Data{9,2}=num2str(round(fout.s2,2));

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
    disp('done');
end

%% Process Fits
[cface1,cedge1] = ixoncolororder(1);



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
errorbar(xvals,phis(:,1)/(2*pi),phis(:,2)/(2*pi),'marker','o',...
    'MarkerFacecolor',cface1,'markeredgecolor',cedge1,'linestyle','none',...
    'linewidth',1.5,'color',cedge1);
% errorbar(xvals,mod(phis(:,1)/(2*pi),1),phis(:,2)/(2*pi),'marker','o',...
    % 'MarkerFacecolor',cface1,'markeredgecolor',cedge1,'linestyle','none',...
    % 'linewidth',1.5,'color',cedge1);
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
ylabel('phase (2\pi)');
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
    
    
    
    m0=(max(yp)-min(yp))./range(xp);
    
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

function guess_out=calculateInitialGuess(x,y,z,opts)    
    thetaVec=linspace(opts.Theta(1),opts.Theta(2),90);    
    
    % remove negative data points
    z(z<0) = 0;    
        
    % Calculate center
    Zx=sum(z,1)'/sum(sum(z));xC=sum(Zx.*x);
    Zy=sum(z,2)/sum(sum(z));yC=sum(Zy.*y);  
    
    % Resize the image (makes it quicker)
    sc=0.2;
    z=imresize(z,sc);x=imresize(x,sc);y=imresize(y,sc);  

    % Compute sum contrasts at many rotation angles
    for jj=1:length(thetaVec)        
        Zrot=imrotate(z,thetaVec(jj));        
        Zsum=sum(Zrot,1);
        Zsum=smooth(Zsum,5);
        CC(jj)=sum(abs(diff(Zsum)).^2);    
    end
    
    % Find the rotation angle which maximizes the contrast
    [~,ind]=max(CC);    
    theta=thetaVec(ind);    

    % Find the cloud gaussian radius
    % Caclulate the second moment orthongal to the fringes
    Zrot=imrotate(z,theta,'crop');
    Zsum1=sum(Zrot,1)/sum(sum(Zrot));
    Zsum2=sum(Zrot,2)/sum(sum(Zrot));   
    yrC=sum(y.*Zsum2);
    sPerp=sqrt(sum(((y-yrC).^2.*Zsum2)));
    
    
    % Assume gaussian radius along fringes is similar
    sParallel = sPerp;
    
    % Find the guess wavelength
    % On the rotated data find the separation of local maxima
    ZsumSmooth=smooth(Zsum1,10);
    [yA,P]=islocalmax(ZsumSmooth,'MinSeparation',50*sc,...
        'MaxNumExtrema',4,'MinProminence',(max(ZsumSmooth)-min(ZsumSmooth))*0.05);
    xA=diff(x(yA));
    L=mean(xA);
    
    if isnan(L) || isinf(L) || isinf(-L) || L<20
       L=80; 
    end
    
    % Find the initial phase of the data
    % Compute correlations of the image with a plane wave at the guess
    % angle
    phiVec=linspace(0,2*pi,50);
    [xx,yy]=meshgrid(x,y);
    Sphi=zeros(length(phiVec),1);
    for nn=1:length(phiVec)        
        plane_wave=sin(2*pi/L*(cosd(theta)*xx+sind(theta)*yy)+phiVec(nn));
        Sphi(nn)=sum(sum(plane_wave.*imgaussfilt(z,1)));     
    end
    
    [~,ind]=max(Sphi);
    phi=phiVec(ind);

    % Guess the modulation strength
    B=max(P)/(max(ZsumSmooth)-min(ZsumSmooth));
    B=0.5;
    
    % Guess the gaussian envelope
    A=sum(sum(z))/(2*pi*sPerp);
    A=max(max(z))/(2*(1+B));
    
    % Construct the initial guess
    guess_out=[A,xC,yC,sPerp,sParallel,B,theta,L,phi];   

end
