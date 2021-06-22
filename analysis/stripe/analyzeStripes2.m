function [hF,outdata]=analyzeStripes2(ixondata,xVar,opts)

global ixon_imgdir
% maskname=fullfile('ixon_mask.mat');
% ixon_mask=load(maskname);
% ixon_mask=ixon_mask.BW;

doDebug=0;

%% Sort the data by the parameter given
params=[ixondata.Params];
xvals=[params.(xVar)];
times=[params.ExecutionDate];

[xvals,inds]=sort(xvals,'ascend');
ixondata=ixondata(inds);

if isequal(xVar,'ExecutionDate')
    xvals=xvals-min(xvals);
end

if nargin==2
   opts.theta=57;
   opts.rotrange=[220 300];
   opts.FitType='Sine';
%    opts.FitType='SmoothSquare';
   opts.LowThreshold=0.2;
   opts.L0=140;
   opts.phi0=pi/2;
   opts.saveAnimation=0;
   opts.B0=0.4;
end

%% Fitting Function
% The fitting function is 2D gaussian who is modulated by a sine wave at
% an angle

gauss2dSine=@(A,xC,yC,sG,B,theta,L,phi,xx,yy) A*...
        exp(-((xx-xC).^2+(yy-yC).^2)/(2*sG^2)).*...
        (1+B*sin(2*pi/L*(cosd(theta)*xx+sind(theta)*yy)+phi));    

myfit=fittype(@(A,xC,yC,sG,B,theta,L,phi,xx,yy) gauss2dSine(A,xC,yC,sG,B,theta,L,phi,xx,yy),...
    'independent',{'xx','yy'},'coefficients',{'A','xC','yC','sG','B','theta',...
    'L','phi'});

opt=fitoptions(myfit);

% Fit tolereances
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
    clf
    thetaVec=linspace(-90,90,100);    
    
    % Get raw data
    Z2=ixondata(kk).Z;
    Z2(Z2<0)=0;
    x2=1:size(Z2,2);x2=x2';
    y2=1:size(Z2,1);y2=y2';    

    % Resize
    sc=0.3;
    Z2=imresize(Z2,sc);
    x2=imresize(x2,sc);
    y2=imresize(y2,sc);
    
    % Fitler
%     Z=imgaussfilt(Z,1);    
    
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
        Zrot=imrotate(Z2,thetaVec(jj),'crop');        
        Zsum=sum(Zrot,1);
        CC(jj)=range(Zsum);
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
    
    if doDebug
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
        waitforbuttonpress
    end
end
t2=now;
disp(['done (' num2str(round((t2-t1)*24*60*60),2) ' s)']);
%% Fit it

fRs={};
hF=figure;
hF.Color='w';
hF.Position(3:4)=[1000 400];
clf
for kk=1:length(ixondata)
    fprintf(['Fitting ' num2str(kk) ' of ' num2str(length(ixondata)) ' ... ']);
    pG=[AG(kk),xCG(kk),yCG(kk),sG(kk),BG(kk),thetaG(kk),LG(kk),phiG(kk)];
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
    
    % Rescale the data for fitting
    sc=0.3;
    Zsc=imresize(Z,sc);
    xsc=imresize(x,sc);
    ysc=imresize(y,sc);
    [xxsc,yysc]=meshgrid(xsc,ysc);
    xxsc(Zsc<10)=[];
    yysc(Zsc<10)=[];
    Zsc(Zsc<10)=[];

    % Create live update figure for fitting    
    figure(hF);
    clf
    colormap(purplemap);

    subplot(231);    
    imagesc(x,y,Z);
    set(gca,'ydir','normal');
    axis equal tight
    colorbar
    
    % Plot the initial guess
    subplot(232)
    imagesc(Zg)    
    set(gca,'ydir','normal');
    axis equal tight
    colorbar
    
    % Perform the fit
    [fout,gof,output]=fit([xxsc(:) yysc(:)],Zsc(:),myfit,opt);
    fRs{kk}=fout;
    % Evalulate the fit
    Zf=feval(fout,xx,yy);
    
    subplot(232)
    imagesc(Zf)    
    set(gca,'ydir','normal');
    colorbar
    
    subplot(233)
    imagesc(Zf-Z);
    set(gca,'ydir','normal');
    colorbar

    subplot(234)
    plot(x,sum(imrotate(Zf,fout.theta,'crop'),1),'r-')
    hold on
    plot(x,sum(imrotate(Z,fout.theta,'crop'),1),'k-'); 
    set(gca,'ydir','normal');
    xlabel('rotated position 1');
    xlim([min(x) max(x)]);
    
    subplot(235)
    plot(y,sum(imrotate(Zf,fout.theta,'crop'),2),'r-')
    hold on
    plot(y,sum(imrotate(Z,fout.theta,'crop'),2),'k-'); 
    set(gca,'ydir','normal');
    xlabel('rotated position 2');
    xlim([min(y) max(y)]);
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


   % Get errors (average +-)
    Ls(kk,2)=(ci(2,7)-ci(1,7))/2;   
    thetas(kk,2)=(ci(2,6)-ci(1,6))/2;    
    phis(kk,2)=(ci(2,8)-ci(1,8))/2;       

end

hF2=figure;
hF2.Position=[100 100 700 250];
hF2.Color='w';
clf;

% Wavelength
subplot(131);
errorbar(xvals,Ls(:,1),Ls(:,2),'marker','o',...
    'MarkerFacecolor',cface1,'markeredgecolor',cedge1,'linestyle','none',...
    'linewidth',1.5,'color',cedge1);
xlabel(xVar);
ylabel('wavelength (px)');

% Angle
subplot(132);
errorbar(xvals,thetas(:,1),thetas(:,2),'marker','o',...
    'MarkerFacecolor',cface1,'markeredgecolor',cedge1,'linestyle','none',...
    'linewidth',1.5,'color',cedge1);
xlabel(xVar);
ylabel('angle (deg.)');

% Phase
subplot(133);
errorbar(xvals,phis(:,1),phis(:,2),'marker','o',...
    'MarkerFacecolor',cface1,'markeredgecolor',cedge1,'linestyle','none',...
    'linewidth',1.5,'color',cedge1);
xlabel(xVar);
ylabel('phase (\pi)');

%% Prepare output data
outdata=struct;
outdata.xVar=xVar;
outdata.xvals=xvals;
outdata.Wavelength=Ls;
outdata.Theta=thetas;
outdata.Phi=phis;
end
