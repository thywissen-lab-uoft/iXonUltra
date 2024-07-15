function [digdata] = dig_compute_radial(digdata,opts)

if nargin == 1
   opts = struct;
end

if ~isfield(opts,'BinStep')
    opts.BinStep = 3;    
end
xcbar = round(mean(digdata.Xc_site));
ycbar = round(mean(digdata.Yc_site));


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
    [rVec,charge,charge_std,n]= radial_profile(Zsub,opts.BinStep);
    
    digdata.rCenter(nn,1) = xc;
    digdata.rCenter(nn,2) = yc;
    digdata.rBinStep = opts.BinStep;
    digdata.r(1:length(rVec),nn) = rVec;
    digdata.Zr(1:length(charge),nn) = charge;
    digdata.Zr_std(1:length(charge),nn) = charge_std;    
    digdata.nr(1:length(charge),nn) = n;
end
    % Send empty indeces to nan for radial vector
    digdata.r([digdata.r==0])=nan;
    digdata.nr(isnan([digdata.r])) = nan;
end

function [Tics,Average,dev,n]=radial_profile(data,radial_step)
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
% dev=accumarray(Z_integer(:),data(:),[],@std);
dev=accumarray(Z_integer(:),data(:),[],@(x) std(x,1));
n= accumarray(Z_integer(:),data(:),[],@(x) numel(x));

end
