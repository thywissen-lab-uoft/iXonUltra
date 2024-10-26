function out = bin_StripeFitCircular(n1,n2,Zb,opts)
%bin_StripeFit Fit a binned fluorescence image to a stripe pattern.
%   When taking fluoresence images of a 2D slice of the 3D cloud,
%   application of a transverse magnetic field induces a stripe pattern to
%  form which is a sample of each plane.  This code fits this distribution
%  and also uses the fluoresence per site to describe the relative focusing
%  between each stripe.
%
%   n1   : lattice site vector which describes the columns of Zb
%   n2   : lattice site vector which describes the rows of Zb
%   Zb   : matrix of size n2xn1 which has fluoresence counts per site
%   opts : options structure
%           - SumIndex : which index are the stripes (NOT WORKING)
%           - LGuess :  wavelength guess
%           - FigNum :  figure number to assign output 
%           - Threshold
if nargin~=4
    opts = struct;
end

% Default figure number
if ~isfield(opts,'FigNum')
    opts.FigNum = 901;
end

% Default wavelength guess is to be automatically determined
if ~isfield(opts,'LGuess')
   opts.LGuess = []; 
end

% Thresholds to throw data away and for focusing score
if ~isfield(opts,'Threshold')
    opts.Threshold = [1000 3000];
end

if ~isfield(opts,'FeshbachField')
    opts.FeshbachField = 124;
end

if ~isfield(opts,'Gradient')
    opts.FeshbachField = 622;
end

[nn1,nn2]=meshgrid(n1,n2);


if ~isfield(opts,'doDebug')
    opts.doDebug=1;
end


%% Wavelength Guess Via 2D FFT
tic

Z = Zb;
Z(isnan(Z))=0;
Z(isinf(Z))=0;

Z(Z<300)=0;
zf = fft2(Z,2^8+1,2^8+1);              % 2D FFT
zf = fftshift(zf);                      % Shift so zero at center
f  = 0.5*linspace(-1,1,size(zf,2));     % Frequency Vector
df=f(2)-f(1);                           % Frequecny spacing
zfnorm=abs(zf);                         % Norm of data

% Radial Profile
[Tics,Average]=radial_profile(zfnorm,1);

% Find peak in radial data
[pks,locs,w,p] =findpeaks(Average,'SortStr','descend','Npeaks',2);

% Get the frequency
f_me = Tics(locs)*df;

% Convert to wavelength

if 1/f_me(1)>100
    ind = 2;
else
    ind =1 ;
end


    guess_lambda= 1/f_me(ind);
lambda_peak_val = pks(ind);
%% Find Rotation Angle
% Find angle with maximum value (at the correct rotation angle, everything
% sums up)
thetaVec = linspace(-pi/4,-pi/4+pi,360);
valTheta = zeros(length(thetaVec),1);
for tt=1:length(thetaVec)
    valTheta(tt)=max(sum(imrotate(zfnorm,thetaVec(tt)*180/pi,'crop'),2));
end
[~,ind]=max(valTheta);

% Rotation Angle
guess_theta = thetaVec(ind);

%% Radius of Curvature
% The radius of curvature depends on the radial field gradient and on the
% bias field.
% Feshbach field --> gotten from uwave frequency
% Bias field gotten from labmda + feshbach field
% Gradient --> gotten from independent measurement
% Radius of curvature --> gotten from bias field and radial gradient

guess_R = 514;

%% Circular Wave Phase
% This function creates a phase map
% L = wavelength
% R      = radius of curvature
% theta  = angle
% phi0   = phase offset

n20 = mean(n2);
y0=n20;
x0=0;

phase_func = @(L,R,theta,phi0,xx,yy) ...
    2*pi/L*(R^2+2*R*(cos(theta)*(xx-x0)+sin(theta)*(yy-y0))+((xx-x0).^2+(yy-y0).^2)).^(1/2)-2*pi*R/L+phi0;

%% Phase Guess

phiVec = linspace(0,2*pi,180);
valPhi = zeros(length(phiVec),1);
for tt=1:length(phiVec)
    map = cos(phase_func(guess_lambda,guess_R,guess_theta,phiVec(tt),nn1,nn2));
    valPhi(tt)=sum(Z.*map,'all');
end
[~,ind]=max(valPhi);

guess_phi=phiVec(ind);
toc

%% Construct Output

guess_phase_map = phase_func(guess_lambda,guess_R,guess_theta,guess_phi,...
    nn1,nn2);
guess_phase_func = @(n1,n2)  phase_func(guess_lambda,guess_R,guess_theta,guess_phi,...
    n1,n2);


%% Histogram All Stripes
Zcopy=Zb;
Zcopy(Zcopy==0)=[];
[n,edges,bin]=histcounts(Zcopy,50);
centers = (edges(1:end-1) + edges(2:end))/2;

%% Kmeans each stripe to determine focusing

phase_map_round = round(guess_phase_map/(2*pi));
stripe_nums = unique(phase_map_round);

histData=zeros(length(n),length(stripe_nums));

Cs=zeros(length(stripe_nums),2);

n1c=zeros(length(stripe_nums),1);
n2c=zeros(length(stripe_nums),1);

for kk=1:length(stripe_nums)
    Zsub=Z.*[phase_map_round==stripe_nums(kk)];

    Zsub(isnan(Zsub))=0;
    Zsub(isinf(Zsub))=0;


   n1c(kk)=sum(nn1.*Zsub,'all')/sum(Zsub,'all');
        n2c(kk)=sum(nn2.*Zsub,'all')/sum(Zsub,'all');

    Zsub(Zsub==0)=[];
    [n,edges,bin]=histcounts(Zsub(:),edges);
    histData(:,kk)=n(:);


    if numel(Zsub)>50
        [idx,c,sumD,D]=kmeans(Zsub(:),2); 
        Cs(kk,1)=min(c);
        Cs(kk,2)=max(c); 

    end
end

sz = 100*Cs(:,2)/max(Cs(:,2));
sz = round(sz);
sz(sz==0)=1;
cc = jet(100);
myc = cc(sz,:);


%% Show it
if opts.doDebug
    FigName = 'BinStripeCircular';
    ff=get(groot,'Children');

    fig=[];
    for kk=1:length(ff)
        if isequal(ff(kk).Name, FigName)
            fig = ff(kk);
        end
    end
    if isempty(fig)
        fig = figure;
    end    
    figure(fig);
    fig.Color='w';
    fig.Name=FigName;
    fig.Position=[50 50 1000 500];
    clf
    
    ax1=subplot(3,4,[1 2 5 6],'parent',fig);
    im1=imagesc(n1,n2,Zb);
    hold on
    colormap(ax1,"parula")
    xlabel('site 1');
    ylabel('site 2');
    axis equal tight    
    ax2=axes;
    ax2.Position=ax1.Position;
    im2=imagesc(n1,n2,cos(guess_phase_map));
    im2.AlphaData=0.2;
    ax2.Visible='off';
    colormap(ax2,"jet")
    axis equal tight  
    linkaxes([ax1,ax2]) 
    ax2.Position=ax1.Position;
    ax1.Position=ax2.Position;
    hold on
    scatter(n1c,n2c,sz,myc,'filled')

    ax_fft=subplot(3,4,[3 4 7 8],'parent',fig);
    imagesc(f,f,zfnorm);
    axis equal tight
    xlim(2*[-1 1]/guess_lambda);
    ylim(2*[-1 1]/guess_lambda);
    colormap(ax_fft,'jet')
    xlabel('fx (1/site)');
    ylabel('fy (1/site)')

    ax_radial=subplot(3,4,9,'parent',fig);
    plot(Tics*df,Average,'.-','parent',ax_radial);
    hold on
    p=plot(1/guess_lambda,lambda_peak_val,'ko','markerfacecolor','k');
    ylim([0 lambda_peak_val]*1.5)
    xlabel('radial frequency (1/site)');
    ylabel('radial average');
    legend(p,{['\lambda = ' num2str(round(guess_lambda,3)) ' site']})
    title('wavelength');

    ax_theta=subplot(3,4,10,'parent',fig);
    plot(thetaVec,valTheta,'.-','parent',ax_theta);
    hold on
    xlabel('rotation angle \theta (rad.)');
    ylabel('sum correlation');
    hold on
    p=plot(guess_theta,valTheta(find(guess_theta==thetaVec,1)),'ko','markerfacecolor','k');
    legend(p,{['\theta = ' num2str(round(guess_theta,3)) ' rad.']})
    title('rotation angle');


    ax_phi=subplot(3,4,11,'parent',fig);
    plot(phiVec,valPhi,'.-','parent',ax_phi);
    hold on
    xlabel('phase angle \phi (rad.)');
    ylabel('sum correlation');
    p=plot(guess_phi,valPhi(find(guess_phi==phiVec,1)),'ko','markerfacecolor','k');
    title('phase');
    legend(p,{['\phi = ' num2str(round(guess_phi,3)) ' rad.']})



   
end

%% output

out = struct;
out.Lambda = guess_lambda;
out.Theta = guess_theta;
out.Radius = guess_R;
out.Phi     = guess_phi;
out.PhaseFunc = @(n1,n2) phase_func(guess_lambda,guess_R,guess_theta,guess_phi,...
    n1,n2);


end



function [Tics,Average]=radial_profile(data,radial_step)
    %main axii cpecified:
    x=(1:size(data,2))-size(data,2)/2;
    y=(1:size(data,1))-size(data,1)/2;
    % coordinate grid:
    [X,Y]=meshgrid(x,y);
    % creating circular layers
    Z_integer=round(abs(X+1i*Y)/radial_step)+1;
    % % illustrating the principle:
    % % figure;imagesc(Z_integer.*data)
    % very fast MatLab calculations:
    Tics=accumarray(Z_integer(:),abs(X(:)+1i*Y(:)),[],@mean);
    Average=accumarray(Z_integer(:),data(:),[],@mean);
end



