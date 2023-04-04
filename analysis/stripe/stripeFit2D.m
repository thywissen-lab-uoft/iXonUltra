function stripeFit2D(z,opt)
%STRIPEFIT2D Summary of this function goes here
%   Detailed explanation goes here

if nargin == 1
   opt = struct;
   opt.FitType = 'sinesquare2';
   opt.doDebug = 0;
   opt.fig = figure;
end

    x = 1:size(z,2);
    y = 1:size(z,1);
    


        % Create a meshgrid
    [xx,yy]=meshgrid(x,y);
%% Do FFT to get initial guess
a = stripeFFT(z,opt);
%% Select fit function
% Define the fitting function. It is a 2D gaussian who is modulated by a
% sine wave at an angle

switch opt.FitType
    case 'sine'

        gauss2dSine=@(A,xC,yC,sG,B,theta,L,phi,xx,yy) A*...
                exp(-((xx-xC).^2+(yy-yC).^2)/(2*sG^2)).*...
                (1+B*cos(2*pi/L*(cos(theta+pi/2)*xx+sin(theta+pi/2)*yy)+phi));    

        myfit=fittype(@(A,xC,yC,sG,B,theta,L,phi,xx,yy) ...
            gauss2dSine(A,xC,yC,sG,B,theta,L,phi,xx,yy),...
            'independent',{'xx','yy'},...
            'coefficients',{'A','xC','yC','sG','B','theta','L','phi'});

        opt=fitoptions(myfit);

        % Define fit options
%         opt.Robust='bisquare';
%         opt.MaxIter=1000;
%         opt.MaxFunEvals=1000;
%         opt.TolFun=1E-8;
%         opt.Robust='bisquare';
%           opt.Robust = 'lar';
% opt.Robust  = 'on';
%         opt.Display='iter';
        
        % Find phase guess (ugh, find a better/quicker way do this
%         phiVec = linspace(0,2*pi,30);
%         c = zeros(length(phiVec),1);
%         
%     % Create the initial guess
%         for kk=1:length(phiVec)            
%             pw = cos(2*pi/a.lambda*(cos(a.theta+pi/2)*xx+sin(a.theta+pi/2)*yy)+phiVec(kk));
%             c(kk) = sum(sum((pw.*z)));
%         end
%         [~,ind]=max(c);
%         
%         a.phi = phiVec(ind);    
        opt.StartPoint = [a.A a.xc a.yc sqrt(a.s1*a.s2) 0.5 a.theta a.lambda a.phi];   
%         opt.Lower = [0 a.xc-5 a.yc-5 sqrt(a.s1*a.s2)*.5 0.3 a.theta-.1 a.lambda*.8 a.phi-1.2*pi];
%         opt.Upper = [a.A*1.5 a.xc+5 a.yc+5 sqrt(a.s1*a.s2)*1.5 0.7 a.theta+.1 a.lambda*1.2 a.phi+1.2*pi];

        Zg=gauss2dSine(a.A,a.xc,a.yc,sqrt(a.s1*a.s2),0.5,a.theta,a.lambda,a.phi,xx,yy); 
        
        case 'sine2'
%%
        gauss2dSine=@(A,theta,L,phi,xx,yy) A*...
                exp(-((xx-a.xc).^2+(yy-a.yc).^2)/(2*a.s1*a.s2)).*...
                0.5.*(1+cos(2*pi/L*(cos(theta+pi/2)*xx+sin(theta+pi/2)*yy)+phi));    

        myfit=fittype(@(A,theta,L,phi,xx,yy) ...
            gauss2dSine(A,theta,L,phi,xx,yy),...
            'independent',{'xx','yy'},...
            'coefficients',{'A','theta','L','phi'});

        opt=fitoptions(myfit);

        % Define fit options
%         opt.Robust='bisquare';
%         opt.MaxIter=1000;
%         opt.MaxFunEvals=1000;
%         opt.TolFun=1E-3; 
%         opt.Robust='lar';
%         opt.Display='iter';
        
        % Find phase guess (ugh, find a better/quicker way do this
%         phiVec = linspace(0,2*pi,30);
%         c = zeros(length(phiVec),1);
%         
%     % Create the initial guess
%         for kk=1:length(phiVec)            
%             pw = cos(2*pi/a.lambda*(cos(a.theta+pi/2)*xx+sin(a.theta+pi/2)*yy)+phiVec(kk));
%             c(kk) = sum(sum((pw.*z)));
%         end
%         [~,ind]=max(c);
%         
%         a.phi = phiVec(ind);        
        
        opt.StartPoint = [a.A a.theta a.lambda a.phi];   
%         opt.Lower = [0  a.theta-.1 a.lambda*.8 a.phi-1.2*pi];
%         opt.Upper = [a.A*1.5  a.theta+.1 a.lambda*1.2 a.phi+1.2*pi];

        Zg=gauss2dSine(a.A,a.theta,a.lambda,a.phi,xx,yy); 
        
        
    
    case 'sinesquare'
        %%
        foo=@(A,xC,yC,sG,eta,theta,L,phi,xx,yy) A*...
                exp(-((xx-xC).^2+(yy-yC).^2)/(2*sG^2)).*...
                0.5.*(1+cos(2*pi/L*(cos(theta+pi/2)*xx+sin(theta+pi/2)*yy)+phi).*...
                1./sqrt(eta^2+cos(2*pi/L*(cos(theta+pi/2)*xx+sin(theta+pi/2)*yy)+phi).^2));

        myfit=fittype(@(A,xC,yC,sG,eta,theta,L,phi,xx,yy) ...
            foo(A,xC,yC,sG,eta,theta,L,phi,xx,yy),...
            'independent',{'xx','yy'},...
            'coefficients',{'A','xC','yC','sG','eta','theta','L','phi'});

        opt=fitoptions(myfit);
%         opt.Robust = 'bisquare';
        
        opt.StartPoint = [a.A a.xc a.yc 1.1*sqrt(a.s1*a.s2) 0.5 a.theta a.lambda a.phi];   
  
        Zg=foo(a.A,a.xc,a.yc,sqrt(a.s1*a.s2),0.5,a.theta,a.lambda,a.phi,xx,yy); 
        
         case 'sinesquare2'
        %%
        foo=@(A,sG,eta,theta,L,phi,xx,yy) A*...
                exp(-((xx-a.xc).^2+(yy-a.yc).^2)/(2*sG^2)).*...
                0.5.*(1+cos(2*pi/L*(cos(theta+pi/2)*xx+sin(theta+pi/2)*yy)+phi).*...
                1./sqrt(eta^2+cos(2*pi/L*(cos(theta+pi/2)*xx+sin(theta+pi/2)*yy)+phi).^2));

        myfit=fittype(@(A,sG,eta,theta,L,phi,xx,yy) ...
            foo(A,sG,eta,theta,L,phi,xx,yy),...
            'independent',{'xx','yy'},...
            'coefficients',{'A','sG','eta','theta','L','phi'});
        opt=fitoptions(myfit);        
        opt.StartPoint = [a.A 1.1*sqrt(a.s1*a.s2) 1 a.theta a.lambda a.phi];    
        Zg=foo(a.A,sqrt(a.s1*a.s2),1,a.theta,a.lambda,a.phi,xx,yy); 
    end
    
        
%% Perform Fit

    figure(fig);
    clf
    subplot(321)
    imagesc(x,y,z)
    axis equal tight
    caxis([0 500]);
    colorbar
    set(gca,'YDir','normal')
    
    
    subplot(323)
    imagesc(x,y,Zg)
     axis equal tight
         set(gca,'YDir','normal')

     
         subplot(324)
         
         
    imagesc(x,y,z-Zg)
     axis equal tight
     colorbar
     caxis([-200 200]);
         set(gca,'YDir','normal')

    % Rescale the data for fitting
    sc=.2;
    Zsc=imresize(z,sc);
    xsc=imresize(x,sc);
    ysc=imresize(y,sc);
    [xxsc,yysc]=meshgrid(xsc,ysc);


    % Perform the fit
    tic
    [fout,gof,output]=fit([xxsc(:) yysc(:)],Zsc(:),myfit,opt);
    toc
    
    disp(a)
    disp(fout)
    disp(gof)

    subplot(325)
    Zf = feval(fout,xx,yy);
    imagesc(Zf);
    axis equal tight
        set(gca,'YDir','normal')

    
    subplot(326)
    imagesc(z-Zf);
    axis equal tight
    colorbar
    caxis([-200 200]);
    set(gca,'YDir','normal')


end

