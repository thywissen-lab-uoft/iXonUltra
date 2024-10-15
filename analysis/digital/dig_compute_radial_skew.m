function [digdata] = dig_compute_radial_skew(digdata,opts)
% dig_compute_radial_skew.m%
%
% Author : C Fujiwara
%
% Despite our best efforts, the digitized data will have some radial skew
% due to imbalances in the chemical potential. (ie. the confining laser
% beams).  We would still to understand the insitu density profile as a
% fuction of the local chemical potential (ie. local density
% approximation).  The only difference is that the equipotentials are now
% ellipses rather than circles.
%


if nargin == 1
   opts = struct;
end

if ~isfield(opts,'BinStep')
    opts.BinStep = 3;    
end
xcbar = round(mean(digdata.Xc_site));
ycbar = round(mean(digdata.Yc_site));

%% Compute Principal Component on average image (for statistics)


R1=[];
R2=[];

for n=1:length(digdata.FileNames)
  R1=[R1 digdata.RatomSite{n}(1,:)];
  R2=[R2 digdata.RatomSite{n}(2,:)];    
end

R=[R1' R2'];
[coeff,score,latent] = pca(R);

% Coefficients are the vectors
v1 = coeff(:,1);    % Semi major axis
v2 = coeff(:,2);    % Semi minor axis

% Eigenvalue is the variance of the vector.  Ordered from most to least
lambda1 = latent(1);
lambda2 = latent(2);

% Calculate ratio of axis
theta = atan2(v1(2),v1(1));    % angle of major axis
ratio = lambda1/lambda2;        % major/minor axis length


%%


    
    digdata=rmfield(digdata,'radial_skew_theta');
    digdata=rmfield(digdata,'radial_skew_ratio');
    digdata=rmfield(digdata,'radial_skew_r');
    digdata=rmfield(digdata,'radial_skew_Zr');
    digdata=rmfield(digdata,'radial_skew_Zr_std');
    digdata=rmfield(digdata,'radial_skew_nr');

for nn=1:size(digdata.Zdig,3)
    % center of this image
    xc = round(digdata.Xc_site(nn));yc = round(digdata.Yc_site(nn));    
    
    if opts.useAverageCenter
       xc = xcbar;
       yc = ycbar;
    end
    % Lattice site vectors
    n1 = digdata.n1;n2 = digdata.n2;        
    nlimits = [min(n1) max(n1) min(n2) max(n2)]; 
    Z = digdata.Zdig(:,:,nn);
    
    %% Compute Square ROI around center of cloud for radial computing    
    L = min(abs(nlimits - [xc xc yc yc]));
    L = L-2;    
    r = [xc xc yc yc]+[-1 1 -1 1]*L;
    
    % Indeces of bounds
    ii = [find(n1 == r(1),1) find(n1 == r(2),1) find(n2 == r(3),1) find(n2 == r(4),1)];
    try
    Zsub = Z(ii(3):ii(4),ii(1):ii(2));
    catch ME
        keyboard
    end
    % Compute radial profile
    [rVec,charge,charge_std,n]= radial_profile_skew(Zsub,opts.BinStep,theta,ratio);
    
    digdata.radial_skew_theta = theta;
    digdata.radial_skew_ratio = ratio;
    digdata.radial_skew_rBinStep = opts.BinStep;
    digdata.radial_skew_r(1:length(rVec),nn) = rVec;
    digdata.radial_skew_Zr(1:length(charge),nn) = charge;
    digdata.radial_skew_Zr_std(1:length(charge),nn) = charge_std;    
    digdata.radial_skew_nr(1:length(charge),nn) = n;
    
end


    % Send empty indeces to nan for radial vector
%     digdata.radial_skew_r([digdata.r==0])=nan;
%     digdata.radial_skew_nr(isnan([digdata.r])) = nan;
end

function [Tics,Average,dev,n]=radial_profile_skew(data,radial_step,theta,ratio)
%main axii cpecified:
x=(1:size(data,2))-size(data,2)/2;
y=(1:size(data,1))-size(data,1)/2;
% coordinate grid:
[X,Y]=meshgrid(x,y);
% creating circular layers

% Map of ellipse where the value is the geometric mean of the major/minor
% axis
geo_mean_map=sqrt((X*cos(theta)+Y*sin(theta)).^2/ratio+ratio*(X*sin(theta)-Y*cos(theta)).^2);


Z_integer = round(geo_mean_map/radial_step)+1;



% Z_integer=round(abs(X+1i*Y)/radial_step)+1;
% % illustrating the principle:
% % figure;imagesc(Z_integer.*data)
% very fast MatLab calculations:
% Tics=accumarray(Z_integer(:),abs(X(:)+1i*Y(:)),[],@mean);
Tics=accumarray(Z_integer(:),geo_mean_map(:),[],@mean);



Average=accumarray(Z_integer(:),data(:),[],@mean);
% dev=accumarray(Z_integer(:),data(:),[],@std);
dev=accumarray(Z_integer(:),data(:),[],@(x) std(x,1));
n= accumarray(Z_integer(:),data(:),[],@(x) numel(x));


smax = 100;
imax=find(Tics>smax,1);
Tics=Tics(1:imax);
Average = Average(1:imax);
dev = dev(1:imax);
n = n(1:imax);



end
