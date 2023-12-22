
%%% Fits stripes to the full image digitization and assigns ROIs


%% Define the data structures

qgmdata_stripes = struct;

%% Define fit function with reasonable bounds

% gauss_square = @(xC,yC,s,L,d,phi,xx,yy) ...
%     (exp(-(...
%     (xx-xC).^2/(2*s^2) + ...
%     (yy-yC).^2/(2*s^2)))).*(square(yy*2*pi/L+phi,d)+1)/2;
% 
% myfit=fittype(@(xC,yC,s,L,d,phi,xx,yy) ...
%     gauss_square(xC,yC,s,L,d,phi,xx,yy),...
%     'independent',{'xx','yy'},...
%     'coefficients',{'xC','yC','s','L','d','phi'});


% gauss_sine = @(xC,yC,L,d,phi,xx,yy) ...
%     (exp(-(...
%     (xx-xC).^2/(2*25^2) + ...
%     (yy-yC).^2/(2*25^2)))).*(sin(yy*2*pi/L+phi)+1)/2;
% 
% myfit=fittype(@(xC,yC,L,d,phi,xx,yy) ...
%     gauss_sine(xC,yC,L,d,phi,xx,yy),...
%     'independent',{'xx','yy'},...
%     'coefficients',{'xC','yC','L','d','phi'});

erf_sine = @(xC,yC,s,r,L,phi,xx,yy) ...
    (1+erf((35^2-((xC-xx).^2/2 + ...
     (yC-yy).^2))/35^2))/2.*(sin(yy*2*pi/L+phi)+1)/2;

myfit=fittype(@(xC,yC,s,r,L,phi,xx,yy) ...
    erf_sine(xC,yC,s,r,L,phi,xx,yy),...
    'independent',{'xx','yy'},...
    'coefficients',{'xC','yC','s','r','L','phi'});


% heavi_square = @(xC,yC,s,L,d,phi,xx,yy) ...
%     heaviside(...
%     s^2 - ((xx-xC).^2 + ...
%     (yy-yC).^2 )).*(square(yy*2*pi/L+phi,d)+1)/2;
% 
% myfit=fittype(@(xC,yC,s,L,d,phi,xx,yy) ...
%     heavi_square(xC,yC,s,L,d,phi,xx,yy),...
%     'independent',{'xx','yy'},...
%     'coefficients',{'xC','yC','s','L','d','phi'});

opts=fitoptions(myfit); 


%% Initial guess

% % xC,yC,L,phi
% opts.StartPoint = [105 85 26 pi];
% opts.Upper = [180 180 180 inf ];
% opts.Lower = [0 0 0 -inf];

% xC,yC,s,r,L,phi
opts.StartPoint = [105 85 35 35 26 pi];
opts.Upper = [180 180 50 180 180 inf ];
opts.Lower = [0 0 1 10 0 -inf];



for ii=1:length(ixondata)

    %%  Perform the Fit

    LD = ixondata(ii).LatticeDig(1);

    xx = LD.n1;
    yy = LD.n2;
    Z  = double(LD.Zdig);
    Natoms_tot = LD.Natoms;

    [xxx,yyy]=meshgrid(xx,yy);

%     % Display the guess
%     fprintf('(xC,yC,L,phi):');
%     s = ['(' ...
%         num2str(round(opts.StartPoint(1))) ',' ...
%         num2str(round(opts.StartPoint(2))) ',' ...
%         num2str(round(opts.StartPoint(3))) ',' ...
%         num2str(round(opts.StartPoint(4))) ')'];
%     fprintf(s);
%     fprintf('->');
% 
%     % Do the fit + time it
%     tic;
%     [fout,gof,output]=fit([xxx(:) yyy(:)],Z(:),myfit,opts);    
%     t=toc;    
%     
%     % Display the fit resut and time to do it
%     s = ['(' ...
%         num2str(round(fout.xC)) ',' ...
%         num2str(round(fout.yC)) ',' ...
%         num2str(round(fout.L,2)) ',' ...
%         num2str(round(fout.phi,2)) ') '];
%     fprintf(s)
%     disp(['(' num2str(t,2) 's) ']);

    % Display the guess
    fprintf('(xC,yC,s,r,L,phi):');
    s = ['(' ...
        num2str(round(opts.StartPoint(1))) ',' ...
        num2str(round(opts.StartPoint(2))) ',' ...
        num2str(round(opts.StartPoint(3))) ',' ...
        num2str(round(opts.StartPoint(4))) ',' ...
        num2str(round(opts.StartPoint(5))) ',' ...
        num2str(round(opts.StartPoint(6))) ')'];
    fprintf(s);
    fprintf('->');

    % Do the fit + time it
    tic;
    [fout,gof,output]=fit([xxx(:) yyy(:)],Z(:),myfit,opts);    
    t=toc;    
    
    % Display the fit resut and time to do it
    s = ['(' ...
        num2str(round(fout.xC)) ',' ...
        num2str(round(fout.yC)) ',' ...
        num2str(round(fout.s)) ',' ...
        num2str(round(fout.r)) ',' ...
        num2str(round(fout.L,2)) ',' ...
        num2str(round(fout.phi,2)) ') '];
    fprintf(s)
    disp(['(' num2str(t,2) 's) ']);

    % Evalulate the fit
    Zfit=feval(fout,xxx,yyy);
    
    err = confint(fout);

%     stipe_fit = struct;
% 
%     stripe_fit.fit = fout;
%     stripe_fit.xC = fout.xC;
%     stripe_fit.yC = fout.yC;
%     stripe_fit.L = fout.L;
%     stripe_fit.phi = fout.phi;
% 
%     stripe_fit.xC = abs(err(1)-fout.xC);
%     stripe_fit.yC = abs(err(3)-fout.yC);
%     stripe_fit.L = abs(err(5)-fout.L);
%     stripe_fit.phi = abs(err(7)-fout.phi);
% 
%     qgmdata_stripes(ii).stripe_fit = stripe_fit;

    %% Assign stripe ROIs + distribute digitization info

    % take cut along xC
    ixC = find(xx==round(fout.xC));

    % find location of fit where phase is pi/2
    diff = -abs(wrapTo2Pi(yy*2*pi/fout.L + fout.phi) - pi/2);
    [p,ip] = findpeaks(diff);

    Zpeaks = Zfit(ip,ixC);
    
    %choose the significant stripes
    threshold = 0.1;
    ip = ip(Zpeaks>threshold);



    
    %Define size of ROIs around center points
    xS = 70;

    for jj=1:length(ip)

        LatticeDig = struct;

        if jj>1 && ceil(yy(ip(jj))-fout.L/2)==floor(yy(ip(jj-1))+fout.L/2)

            %make sure ROIs don't overlap so atoms aren't double counted
            ROI = [xx(ixC)-xS xx(ixC)+xS ceil(yy(ip(jj))-fout.L/2)+1 floor(yy(ip(jj))+fout.L/2)];
        else
            ROI = [xx(ixC)-xS xx(ixC)+xS ceil(yy(ip(jj))-fout.L/2) floor(yy(ip(jj))+fout.L/2)];
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

    

    %% Assess the fit

    LatticeDig = qgmdata_stripes(ii).LatticeDig; 
    NStripes = sum([LatticeDig(:).Natoms]);

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
    
    % Residue plot
    ax3=subplot(1,2,2);    
    hImg_err=imagesc(xx,yy,Zfit-Z);
    hold on;
    plot(fout.xC,fout.yC,'rx','MarkerSize',18)
    plot(xx(ixC)*ones(length(ip),1),yy(ip),'ok')

    for i=1:length(ip)
        plot(LatticeDig(i).ROI(1:2),LatticeDig(i).ROI(3)*ones(2,1),'r')
        plot(LatticeDig(i).ROI(1:2),LatticeDig(i).ROI(4)*ones(2,1),'c')
    end
    set(gca,'ydir','normal');
    axis equal tight
    colorbar;


keyboard

    
end




