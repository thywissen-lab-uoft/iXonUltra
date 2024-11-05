function out = StripeCircle(X,Y,Z,opts)


[XX,YY]=meshgrid(X,Y);

if nargin==3
    opts=struct;
end

if ~isfield(opts,'doDebug')
    opts.doDebug=1;
end

if ~isfield(opts,'Name')
    opts.Name=[];
end

disp('pixel domain stripe analysis');


% 16 um at a magnification of 
a_site = 2.68;

%% Get Data And FFT

Z(isnan(Z))=0;
Z(isinf(Z))=0;

zf = fft2(Z,2^12+1,2^12+1);              % 2D FFT
zf = fftshift(zf);                      % Shift so zero at center
f  = 0.5*linspace(-1,1,size(zf,2));     % Frequency Vector
df=f(2)-f(1);                           % Frequecny spacing
zfnorm=abs(zf);                         % Norm of data
%% Find Rotation Angle

% Mininum lambda
lambda_min = 10*a_site;

px_max = ceil((1/lambda_min)/df);
px_center = (size(zfnorm,1)+1)/2;
fmax = f(px_center+px_max);

px_lims = (px_center-px_max):(px_center+px_max);
f_sub = f((px_center-px_max):(px_center+px_max));
[ffx_sub,ffy_sub]=meshgrid(f_sub,f_sub);
fmat=(ffx_sub.^2+ffy_sub.^2).^(1/2);
F = fmat<=fmax;

zfnorm_sub = zfnorm(px_lims,px_lims).*F;

% Find angle with maximum value (at the correct rotation angle, everything
% sums up)
thetaVec = linspace(-pi/4,-pi/4+pi,360);
valTheta = zeros(length(thetaVec),1);
tic
for tt=1:length(thetaVec)
    % valTheta(tt)=max(sum(imrotate(zfnorm,thetaVec(tt)*180/pi,'crop'),2));
    valTheta(tt)=max(sum(imrotate(zfnorm_sub,thetaVec(tt)*180/pi,'crop'),2));
end
tRotate=toc;
[valTheta0,ind]=max(valTheta);

inds=[valTheta>0.9*valTheta0];


% Rotation Angle
% guess_theta = thetaVec(ind);
guess_theta = sum(valTheta(inds).*thetaVec(ind))/sum(valTheta(inds));
guess_theta = median(thetaVec(inds));

%% Wavelength 

zfnorm_rotate = imrotate(zfnorm,180/pi*guess_theta,'crop');

iCenter = (size(zfnorm_rotate,1)+1)*0.5;
zfnorm_rotate_sub = zfnorm_rotate((iCenter-10:iCenter+10),:);

zfnorm_rotate_sub_sum=sum(zfnorm_rotate_sub,1);

lambdaMin = 40; % in pixels
lambdaMax = 120;

fMin = 1/lambdaMax;
fMax = 1/lambdaMin;

[pks,locs,w,p] = findpeaks(zfnorm_rotate_sub_sum,f,'NPeaks',3,'SortStr',...
    'descend','MinPeakDistance',fMin);

[locs_sort,inds]=sort(locs,'ascend');
fMe = 0.5*(locs_sort(3)-locs_sort(1));

pks_sort = pks(inds);

lambda_peak_val = mean(pks_sort(1:2));
guess_lambda = 1/fMe;

% 
% 
% % Radial Profile
% [Tics,Average]=radial_profile(zfnorm,1);
% 
% % Find peak in radial data
% [pks,locs,w,p] =findpeaks(Average,'SortStr','descend','Npeaks',2);
% 
% 
% % Get the frequency
% f_me = Tics(locs)*df;
% 
% % Convert to wavelength
% 
% if 1/f_me(1)>100
%     ind = 2;
% else
%     ind =1 ;
% end
% 
% guess_lambda= 1/f_me(ind);
% lambda_peak_val = pks(ind);

%% Radius of Curvature
% The radius of curvature depends on the radial field gradient and on the
% bias field.
% Feshbach field --> gotten from uwave frequency
% Bias field gotten from labmda + feshbach field
% Gradient --> gotten from independent measurement
% Radius of curvature --> gotten from bias field and radial gradient

guess_R = 514*a_site;

%% Circular Wave Phase
% This function creates a phase map
% L = wavelength
% R      = radius of curvature
% theta  = angle
% phi0   = phase offset

y0=mean(Y);
x0=0;

phase_func = @(L,R,theta,phi0,xx,yy) ...
    2*pi/L*(R^2+2*R*(cos(theta)*(xx-x0)+sin(theta)*(yy-y0))+((xx-x0).^2+(yy-y0).^2)).^(1/2)-2*pi*R/L+phi0;

%% Phase Guess
% Compute best phase using phasors
Cmap = sum(Z.*cos(phase_func(guess_lambda,guess_R,guess_theta,0,XX,YY)),'all');
Smap = sum(Z.*cos(phase_func(guess_lambda,guess_R,guess_theta,pi/2,XX,YY)),'all');
guess_phi=atan2(Smap,Cmap);


%% Construct Output

guess_phase_map = phase_func(guess_lambda,guess_R,guess_theta,guess_phi,...
    XX,YY);
guess_phase_func = @(x,y)  phase_func(guess_lambda,guess_R,guess_theta,guess_phi,...
    x,y);


%% Show it
if opts.doDebug

    if ~isfield(opts,'Parent')
    FigName = 'StripeCircular';
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
    % fig.ToolBar='none';
    % fig.MenuBar='none';
    fig.Position=[5 50 350 300];
    clf(fig);
    else
        fig = opts.Parent;
        for kk=1:length(fig.Children)
            delete(fig.Children(1))
        end
    end
    ca = [1 1 1];
    cb = [0.6 0 .5];
    cAtoms = [linspace(ca(1),cb(1),1000)' ...
        linspace(ca(2),cb(2),1000)' linspace(ca(3),cb(3),1000)'];

    ca = [1 1 1];
    cb = [255,215,0]/255;
    cStripe = [linspace(ca(1),cb(1),1000)' ...
        linspace(ca(2),cb(2),1000)' linspace(ca(3),cb(3),1000)'];

    if ~isempty(opts.Name)
        t=uicontrol('parent',fig,'units','pixels','string',opts.Name,...
            'fontsize',8,'horizontalalignment','left','backgroundcolor','w',...
            'style','text');
        t.Position=[2 2 300 15];
    end

    ax1=subplot(3,1,[1 2],'parent',fig);
    imagesc(X,Y,Z,'parent',ax1);
    colormap(ax1,cAtoms);
    xlabel(ax1,'x');
    ylabel(ax1,'y');
    set(ax1,'YDir','normal');
    caxis(ax1,[0 max(Z,[],'all')*0.2]);
    axis(ax1,'tight');axis(ax1,'equal'); 
    hold(ax1,'on');
    ax2=axes('parent',fig);
    ax2.Position=ax1.Position;
    im2=imagesc(X,Y,cos(guess_phase_map),'parent',ax2);
    im2.AlphaData=0.3;
    ax2.Visible='off';
    colormap(ax2,cStripe)
    axis(ax2,'tight');axis(ax2,'equal'); 
    linkaxes([ax1,ax2]) 
    ax2.Position=ax1.Position;
    ax1.Position=ax2.Position;
    set(ax1,'fontsize',10);
    set(ax2,'fontsize',10);  
    set(ax2,'YDir','normal');
    title(ax1,'image');

    ax_fft=subplot(3,3,7,'parent',fig);
    imagesc(f,f,zfnorm,'parent',ax_fft);
    axis(ax_fft,'equal');
    axis(ax_fft,'tight');
    xlim(ax_fft,2*[-1 1]/guess_lambda);
    ylim(ax_fft,2*[-1 1]/guess_lambda);
    colormap(ax_fft,'jet')
    xlabel(ax_fft,'fx (1/px)');
    ylabel(ax_fft,'fy (1/px)')
    set(ax_fft,'YDir','normal','fontsize',10);
    title(ax_fft,'abs fft');

    % ax_radial=subplot(3,3,8,'parent',fig);
    % plot(Tics*df,Average,'.-','parent',ax_radial);
    % ylim(ax_radial,[0 lambda_peak_val]*1.5)
    % xlabel(ax_radial,'radial frequency (1/px)');
    % ylabel(ax_radial,'radial average');
    % hold(ax_radial,'on');
    % p=plot([1 1]/guess_lambda,[0 1]*lambda_peak_val,'k-','markerfacecolor','k',...
    %     'parent',ax_radial);    
    % legend(p,{['\lambda = ' num2str(round(guess_lambda,3)) ' px']},...
    %     'parent',ax_radial.Parent);
    % xlim(ax_radial,[0 0.2])
    % set(ax_radial,'fontsize',10);
    % title(ax_radial,'wavelength');

    ax_sum=subplot(3,3,8,'parent',fig);
    plot(f,zfnorm_rotate_sub_sum,'.-','parent',ax_sum);
    ylim(ax_sum,[0 lambda_peak_val]*1.5)
    xlabel(ax_sum,'rotated frequency (1/px)');
    ylabel(ax_sum,'partial sum profile average');
    hold(ax_sum,'on');
    p=plot([1 1]/guess_lambda,[0 1]*lambda_peak_val,'k-','markerfacecolor','k',...
        'parent',ax_sum);  
    p=plot(-[1 1]/guess_lambda,[0 1]*lambda_peak_val,'k-','markerfacecolor','k',...
        'parent',ax_sum);    
    legend(p,{['\lambda = ' num2str(round(guess_lambda,3)) ' px']},...
        'parent',ax_sum.Parent);
    xlim(ax_sum,[-2.5 2.5]/guess_lambda)
    set(ax_sum,'fontsize',10);
    title(ax_sum,'wavelength');


    

    ax_theta=subplot(3,3,9,'parent',fig);
    plot(180/pi*thetaVec,valTheta,'.-','parent',ax_theta);
    xlabel(ax_theta,'rotation angle \theta (deg.)');
    ylabel(ax_theta,'sum correlation');
    set(ax_theta,'Ylim',[min(valTheta) max(valTheta)*1.1]);
    hold(ax_theta,'on');
    p=plot(180/pi*guess_theta*[1 1],[min(valTheta) max(valTheta)],'k-','markerfacecolor','k','parent',ax_theta);
    legend(p,{['\theta = ' num2str(round(180/pi*guess_theta,1)) ' deg.']},...
        'location','best','parent',ax_theta.Parent)
    set(ax_theta,'fontsize',10);
    title(ax_theta,'rotation');


end

%% output

out = struct;
out.Lambda = guess_lambda;
out.Theta = guess_theta;
out.Radius = guess_R;
out.Phi     = guess_phi;
out.PhaseFunc = @(x,y) phase_func(guess_lambda,guess_R,guess_theta,guess_phi,...
    x,y);




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

