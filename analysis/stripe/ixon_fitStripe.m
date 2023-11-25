function stripe=ixon_fitStripe(data,opts)

if nargin ==1
    opts = struct;
    % opts.xVar = 'ExecutionDate';
end

if ~isfield(opts,'Theta');
   opts.Theta=[10 190]; 
end

sc = 0.3;

%% Fitting Function
% https://en.wikipedia.org/wiki/Gaussian_function

% Define the fitting function as a 2D rotated gaussian whose principal axes
% are the same as the plane wave modulation

a = @(theta,s1,s2) cosd(theta)^2/(2*s1^2) + sind(theta)^2/(2*s2^2);
b = @(theta,s1,s2) -sind(2*theta)/(4*s1^2) + sind(2*theta)/(4*s2^2);
c = @(theta,s1,s2) sind(theta)^2/(2*s1^2) + cosd(theta)^2/(2*s2^2);
phase_map = @(L,theta,phi,xx,yy) ...
    2*pi/L*(cosd(theta)*xx+sind(theta)*yy)+phi;
rot_gauss = @(A,xC,yC,s1,s2,theta,xx,yy) ...
    A.*exp(-(...
    a(theta,s1,s2)*(xx-xC).^2 + ...
    2*b(theta,s1,s2)*(xx-xC).*(yy-yC)+ ...
    c(theta,s1,s2)*(yy-yC).^2));
gauss2dSineRot = @(A,xC,yC,s1,s2,B,theta,L,phi,xx,yy) ...
    rot_gauss(A,xC,yC,s1,s2,theta,xx,yy).*...
    (1+B*sin(phase_map(L,theta,phi,xx,yy)));
myfit2=fittype(@(A,xC,yC,s1,s2,B,theta,L,phi,xx,yy) ...
    gauss2dSineRot(A,xC,yC,s1,s2,B,theta,L,phi,xx,yy),...
    'independent',{'xx','yy'},...
    'coefficients',{'A','xC','yC','s1','s2','B','theta','L','phi'});
opt2=fitoptions(myfit2);             

%% Grab the Data
      
Z=data.Z;
Z=data.ZNoFilter;
Z=imgaussfilt(Z,1);
x=data.X;x=x';
y =data.Y;y=y'; 

%% Create Initial Guess
    
tic;
G=calculateInitialGuess(x,y,Z,opts);
t=toc;

opt2.StartPoint = G;

% Assign reasonable upper and lower bounds
% A,xC,yC,s1,s2,B,theta,L,phi
opt2.Upper = [G(1)*2 512 512 250 250 1  360 512 inf];
opt2.Lower = [0        0   0   0   0 0 -360   0 -inf];
%% Rescale the data

Zsc=imresize(Z,sc);xsc=imresize(x,sc);ysc=imresize(y,sc);
[xxsc,yysc]=meshgrid(xsc,ysc);

%% Remove signal less than or equal to zero.   

xxsc(Zsc<=0)=[];yysc(Zsc<=0)=[];Zsc(Zsc<=0)=[];        
    
%% Remove Data points outside the ROI

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

%% Perform the Fit

% Display the guess
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

% Display the fit resut and time to do it
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

%% Return the output

stripe = struct;

% Fit output and objects
stripe.Fit = fout;
stripe.gof = gof;
stripe.rsquare = gof.rsquare;
stripe.sse = gof.sse;

stripe.PhaseMapFunc = phase_map;
stripe.EnvelopeFunc = rot_gauss;

ci=confint(fout);

% Fit coefficients and errors
stripe.A = fout.A;
stripe.A_err = (ci(2,1)-ci(1,1))*.5;

stripe.xC = fout.xC;
stripe.xC_err = (ci(2,2)-ci(1,2))*.5;

stripe.yC = fout.yC;
stripe.yC_err = (ci(2,3)-ci(1,3))*.5;

stripe.s1 = fout.s1;
stripe.s1_err = (ci(2,4)-ci(1,4))*.5;

stripe.s2 = fout.s2;
stripe.s2_err = (ci(2,5)-ci(1,5))*.5;


stripe.B = fout.B;
stripe.B_err = (ci(2,6)-ci(1,6))*.5;

stripe.theta = fout.theta;
stripe.theta_err = (ci(2,7)-ci(1,7))*.5;

stripe.L = fout.L;
stripe.L_err = (ci(2,8)-ci(1,8))*.5;

stripe.phi = fout.phi;
stripe.phi_err = (ci(2,9)-ci(1,9))*.5;


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
