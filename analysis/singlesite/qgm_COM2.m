%%% Fits stripes to the full image digitization by summing along the x direction 
% and assigns ROIs


%% Define the data structures

qgmdata_stripes = struct;

%% Define sinusoidal fit function in 1D 

% erf_sine_1D = @(yC,s,A,L,phi,yy) ...
%     (1+erf((s^2 - (yC-yy).^2)/r^2))/2.*(sin(yy*2*pi/L+phi)+1)/2;

gauss_sine_1D = @(yC,s,A,L,phi,m,yy) ...
    A*(exp(-(...
    (yy-yC).^2/(2*s^2)))).*(m*sin(yy*2*pi/L+phi)+1)/2;

myfit_1D=fittype(@(yC,s,A,L,phi,m,yy) ...
    gauss_sine_1D(yC,s,A,L,phi,m,yy),...
    'independent',{'yy'},...
    'coefficients',{'yC','s','A','L','phi','m'});

opts_1D=fitoptions(myfit_1D); 

%% Define 2D fit function

gauss_sine = @(xC,yC,L,sX,sY,phi,xx,yy) ...
    (exp(-(...
    (xx-xC).^2/(2*sX^2) + ...
    (yy-yC).^2/(2*sY^2)))).*(sin(yy*2*pi/L+phi)+1)/2;

erf_sine = @(xC,s,r,yC,L,phi,xx,yy) ...
    (1+erf((s^2-((xC-xx).^2/2 + ...
     (yC-yy).^2))/r^2))/2.*(sin(yy*2*pi/L+phi)+1)/2;

myfit=fittype(@(xC,s,r,yC,L,phi,xx,yy) ...
    erf_sine(xC,s,r,yC,L,phi,xx,yy),...
    'independent',{'xx','yy'},...
    'coefficients',{'xC','s','r','yC','L','phi',});

opts=fitoptions(myfit); 


%% Initial guess

% yC,s,r,L,phi,m
opts_1D.StartPoint = [LD.Yc LD.Ys 0.05 26 pi 1];
opts_1D.Upper = [180 50 1 180 inf 1];
opts_1D.Lower = [0 1  0 0 -inf 0];
% % yC,s,r,L,phi,m,d
% opts_1D.StartPoint = [LD.Yc LD.Ys 0.05 26 pi 1 0.5];
% opts_1D.Upper = [180 50 1 180 inf 1 1 ];
% opts_1D.Lower = [0 1  0 -inf 0 0];

% xC,s,r
opts.StartPoint = [105 35 35 0 0 0 ];
opts.Upper = [180 50 180 0 0 0 ];
opts.Lower = [0 0 10 0 0 0 ];


for ii=1:length(ixondata)

    %% Load the data

    LD = ixondata(ii).LatticeDig(1);

    xx = LD.n1;
    yy = LD.n2;
    Z  = double(LD.Zdig);
    Natoms_tot = LD.Natoms;

    [xxx,yyy]=meshgrid(xx,yy);

    % Sum Z data along X direction to get initial sine fit
    Zy = sum(Z,2)/Natoms_tot;
    Zx = sum(Z,1);

    % Smooth the data with 7 neighboring sites (3 on each side)
    Zy = smooth(Zy,5);

    %% Perform 1D fit
    % We fit the summed counts along the stripe direction to find the
    % wavelength of the stripe pattern

        % Display the guess
    fprintf('1D (yC,s,r,L,phi):');
    s = ['(' ...
        num2str(round(opts_1D.StartPoint(1))) ',' ...
        num2str(round(opts_1D.StartPoint(2))) ',' ...
        num2str(round(opts_1D.StartPoint(3),2)) ',' ...
        num2str(round(opts_1D.StartPoint(4))) ',' ...
        num2str(round(opts_1D.StartPoint(5))) ',' ...
        num2str(round(opts_1D.StartPoint(6),2)) ')'];
    fprintf(s);
    fprintf('->');

    % Do the fit + time it
    tic;
    [fout_1D,gof_1D,output_1D]=fit(yy(:),Zy(:),myfit_1D,opts_1D);    
    t=toc;  

    % Display the fit resut and time to do it
    s = ['(' ...
        num2str(round(fout_1D.yC)) ',' ...
        num2str(round(fout_1D.s)) ',' ...
        num2str(round(fout_1D.A,2)) ',' ...
        num2str(round(fout_1D.L,2)) ',' ...
        num2str(round(fout_1D.phi,2)) ',' ...
        num2str(round(fout_1D.m,2)) ') '];
    fprintf(s)
    disp(['(' num2str(t,2) 's) ']);

    % Evalulate the fit
    Zfit_1D=feval(fout_1D,yy);
    err_1D = confint(fout_1D);
    

    % Plot the fit

    fig_1d = figure(3001);
    plot(yy,Zy,'.')
    hold on;
    plot(yy,Zfit_1D)
    hold off;

    %% Old 2D fit
%     % Update guess of size based off 1D fit
%     opts.StartPoint(2) = 2*fout_1D.s;
% 
%     % Constrain fit variables from 1D fit
%     opts.StartPoint(4) = fout_1D.yC;
%     opts.StartPoint(5) = fout_1D.L;
%     opts.StartPoint(6) = fout_1D.phi;
% 
%     opts.Lower(2) = fout_1D.L;
%     opts.Lower(4) = fout_1D.yC;
%     opts.Lower(5) = fout_1D.L;
%     opts.Lower(6) = fout_1D.phi;
% 
%     opts.Upper(4) = fout_1D.yC;
%     opts.Upper(5) = fout_1D.L;
%     opts.Upper(6) = fout_1D.phi;
% 
%     % Display the guess
%     fprintf('2D (xC,s,r):');
%     s = ['(' ...
%         num2str(round(opts.StartPoint(1))) ',' ...
%         num2str(round(opts.StartPoint(2))) ',' ...
%         num2str(round(opts.StartPoint(3))) ')'];
%     fprintf(s);
%     fprintf('->');
% 
%     % Do the fit + time it
%     tic;
%     [fout,gof,output]=fit([xxx(:) yyy(:)],Z(:),myfit,opts);    
%     t=toc;    
%     
%     % Display the fit result and time to do it
%     s = ['(' ...
%         num2str(round(fout.xC)) ',' ...
%         num2str(round(fout.s)) ',' ...
%         num2str(round(fout.r))  ') '];
%     fprintf(s)
%     disp(['(' num2str(t,2) 's) ']);
% 
%     % Evalulate the fit
%     
%     Zfit=feval(fout,xxx,yyy);
%    
%     
%     
%     err = confint(fout);
    
%% Construct 2D expression

Zfit = gauss_sine(LD.Xc,LD.Yc,fout_1D.L,LD.Xs,LD.Ys,fout_1D.phi,xxx,yyy);

    
    stipe_fit = struct;

%     stripe_fit.fit = fout;
%     stripe_fit.xC = fout.xC;
%     stripe_fit.yC = fout_1D.yC;
%     stripe_fit.L = fout_1D.L;
%     stripe_fit.phi = fout_1D.phi;
% 
%     stripe_fit.xC_err = abs(err(1)-fout.xC);
%     stripe_fit.yC_err = abs(err_1D(1)-fout_1D.yC);
%     stripe_fit.L_err = abs(err_1D(7)-fout_1D.L);
%     stripe_fit.phi_err = abs(err_1D(9)-fout_1D.phi);

    stripe_fit.fit = fout;
    stripe_fit.xC = LD.Xc;
    stripe_fit.yC = fout_1D.yC;
    stripe_fit.L = fout_1D.L;
    stripe_fit.phi = fout_1D.phi;
    stripe_fit.mod_depth = fout_1D.m;

    stripe_fit.yC_err = abs(err_1D(1)-fout_1D.yC);
    stripe_fit.L_err = abs(err_1D(7)-fout_1D.L);
    stripe_fit.phi_err = abs(err_1D(9)-fout_1D.phi);
    stripe_fit.mod_depth_err = abs(err_1D(11)-fout_1D.m);

    qgmdata_stripes(ii).stripe_fit = stripe_fit;

        %% Assign stripe ROIs + distribute digitization info

    % take cut along xC
    ixC = find(xx==round(stripe_fit.xC));

    % find location of fit where phase is pi/2
    diff = -abs(wrapTo2Pi(yy*2*pi/stripe_fit.L + stripe_fit.phi) - pi/2);
    [p,ip] = findpeaks(diff);

    Zpeaks = Zfit(ip,ixC);

%     threshold = 0.05;
%     ip = ip(Zpeaks>threshold);
    
    
    %choose the 4 most significant stripes
    [Ztop, itop] = maxk(Zpeaks,4);
     ip = ip(sort(itop));


    
    %Define size of ROIs around center points
    xS = 70;

    for jj=1:length(ip)

        LatticeDig = struct;

        if jj>1 %&& ceil(yy(ip(jj))-stripe_fit.L/2)==floor(yy(ip(jj-1))+stripe_fit.L/2)

            %make sure ROIs don't overlap so atoms aren't double counted
            ROI = [xx(ixC)-xS xx(ixC)+xS floor(yy(ip(jj-1))+stripe_fit.L/2)+1 floor(yy(ip(jj))+stripe_fit.L/2)];
        else
            ROI = [xx(ixC)-xS xx(ixC)+xS ceil(yy(ip(jj))-stripe_fit.L/2) floor(yy(ip(jj))+stripe_fit.L/2)];
        end

        i1 = find(xx==ROI(1));
        i2 = find(xx==ROI(2));
        i3 = find(yy==ROI(3));
        i4 = find(yy==ROI(4));

        xStripe = xx(i1:i2);
        yStripe = yy(i3:i4);
        ZStripe = Z(i3:i4,i1:i2);

        Natoms = sum(sum(ZStripe));        % Total number of atoms


        zY=sum(ZStripe,2)';zY = zY/sum(zY);
        zX=sum(ZStripe,1); zX = zX/sum(zX);

        % Calculate center of mass
        Xc=sum(zX.*xStripe);
        Yc=sum(zY.*yStripe);          

        % Calculate central second moment/variance and the standard
        % deviation
        X2=sum(zX.*(xStripe-Xc).^2); % x variance
        Xs=sqrt(X2); % standard deviation X
        Y2=sum(zY.*(yStripe-Yc).^2); % x variance
        Ys=sqrt(Y2); % standard deviation Y               

        LatticeDig.ROI  = ROI;
        LatticeDig.Zdig = ZStripe;
        LatticeDig.Natoms = Natoms;
        LatticeDig.n1 = xStripe;
        LatticeDig.n2 = yStripe;
        LatticeDig.Xc = Xc;
        LatticeDig.Yc = Yc;
        LatticeDig.Xs = Xs;
        LatticeDig.Ys = Ys;  
        

        qgmdata_stripes(ii).LatticeDig(jj) = LatticeDig;



    end

    
% 
    %% Assess the fit

    LatticeDig = qgmdata_stripes(ii).LatticeDig; 
    NStripes = sum([LatticeDig(:).Natoms]);

    fprintf('%.f stripes above threshold \n',length(ip))
    fprintf('Total image atom number: %.f \n',Natoms_tot)
    fprintf('Total stripe atom number: %.f \n',NStripes)

    if Natoms_tot<NStripes
        fprintf('%.f atoms double counted (%.2f %%) \n',NStripes-Natoms_tot,(NStripes-Natoms_tot)/Natoms_tot*100)
    elseif Natoms_tot>NStripes
        fprintf('%.f atoms missed (%.2f %%) \n',-NStripes+Natoms_tot,-(NStripes-Natoms_tot)/Natoms_tot*100)
    end

%% Show the data and fit


    % Make the figure
    hF_live=figure(2005);
    hF_live.Color='w';
    hF_live.Position=[100 500 1100 500];
    
    % Image Plot
    ax1=subplot(1,2,1);    
    hImg_raw=imagesc(xx,yy,Z);
    set(gca,'ydir','normal');
    axis equal tight  
    xlim([min(xx) max(xx)])
    ylim([min(yy) max(yy)])
    
    % Residue plot
    ax3=subplot(1,2,2);    
    hImg_err=imagesc(xx,yy,Zfit-Z);
    hold on;
    plot(stripe_fit.xC,stripe_fit.yC,'rx','MarkerSize',18)
    plot(xx(ixC)*ones(length(ip),1),yy(ip),'ok')


    for i=1:length(ip)
        plot(LatticeDig(i).ROI(1:2),LatticeDig(i).ROI(3)*ones(2,1),'r')
        plot(LatticeDig(i).ROI(1:2),LatticeDig(i).ROI(4)*ones(2,1),'c')
    end
    set(gca,'ydir','normal');
    axis equal tight
    colorbar;

    xlim([min(xx) max(xx)])
    ylim([min(yy) max(yy)])

    keyboard
end

% Save the output data
filename=fullfile(ixon_imgdir,'figures','qgmdata_stripes.mat');
save(filename,'qgmdata_stripes');

