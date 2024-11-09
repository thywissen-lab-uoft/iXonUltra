function out=binLattice(x,y,z,opts)
% Author : C. Fujiwara
%
% Given a fluoresence image input (z) with corresponding pixel position
% cooridinates (x,y), this function assigns fluoresnce to each lattice
% site.
%
% The phase is found my optizing the total contrast of the image which is
% gotten from the 

% Matrix which defines lattice basis in units of camera pixels
A = [opts.a1 opts.a2];                

% Image pixel size in units of camera pixels
% dX = x(2)-x(1);
% dY = y(2)-y(1);

if ~isfield(opts,'PixelThreshold')
   opts.PixelThreshold = 100; 
end

z(z<opts.PixelThreshold)=0;

%% Image preparation

% Rescale the image to avoid rounding errors when moving pixels
% scale amplitude to preserve norm

z2 = imresize(z,opts.ScaleFactor,'method','bilinear')/(opts.ScaleFactor^2); 
x2 = linspace(x(1),x(end),size(z2,2));
y2 = linspace(y(1),y(end),size(z2,1));

% Create 2xN vectors of all pixels positions
% tic
[X,Y]=meshgrid(x2,y2);                  % matrix of X and Y
R = [X(:) Y(:)]'; % long part, becaues of transpose?
% toc
% 
if isfield(opts,'PixelThreshold') && isfield(opts,'UsePixelThreshold')
    z2(z2<(opts.PixelThreshold/opts.ScaleFactor^2))=0;    
end

% 1xN vector all counts per pixel
Z = z2(:);                              % counts per re-scaled pixel

%% Convert Pixels to Lattice Sites
% Takes 0.18 seconds at scale factor 10

% Solve for all position in terms of the lattice basis
N0 = inv(A)*R;


p = [opts.p1; opts.p2];             % Phase        
P = repmat(p,[1 length(Z)]);        % Phase vector      

% Remove phase from calculated bin positions and round to the nearest
% lattice site
M = round(N0-P);                    % round pixel to lattice site


%% Find bounds on lattice indeces

% Given the lattice basis and the image bounds, determine some reasonable
% bounds on the extremal lattice indeces.

% Solve for position of each corner in terms of the lattice basis
na = inv(A)*[x(1); y(1)];
nb = inv(A)*[x(1); y(end)];
nc = inv(A)*[x(end); y(1)];
nd = inv(A)*[x(end); y(end)];

% Collect the four corners
Nc = [na nb nc nd];

% Bounds and vector for lattice basis vector 1
n1i = floor(min(Nc(1,:)));
n1f = ceil(max(Nc(1,:)));
n1 = n1i:n1f;

% Bounds and vector for lattice basis vector 2
n2i = floor(min(Nc(2,:)));
n2f = ceil(max(Nc(2,:)));
n2 = n2i:n2f;

%% Final Binning

% Takes 0.3 seconds at x10 scale

% Remove points outside of the desired ROI.
ibad = logical((M(1,:)<n1i) + (M(1,:)>n1f) + (M(2,:)>n2f) + (M(2,:)<n2i));
M(:,ibad) = [];
Z(ibad) = [];

Zbin = zeros(n2f-n2i+1,n1f-n1i+1);
Znum = zeros(n2f-n2i+1,n1f-n1i+1);
for ii=1:size(M,2)    
    m1 = M(1,ii)-n1i+1;
    m2 = M(2,ii)-n2i+1;    
    Znum(m2,m1) = Znum(m2,m1)+1;
    Zbin(m2,m1) = Zbin(m2,m1) + Z(ii);  
end

%% Debugging
doDebug=0;
if doDebug
    [N1,N2]=meshgrid(n1,n2);
    x_site = opts.a1(1).*(N1+opts.p1) + opts.a2(1).*(N2+opts.p2);
    y_site = opts.a1(2).*(N1+opts.p1) + opts.a2(2).*(N2+opts.p2);
    
    hF = figure;
    hF.Color='w';
    
    subplot(121);
    imagesc(x2,y2,z2);   
    hold on
   z3 = zeros(numel(y2),numel(x2));


    plot(x_site(:),y_site(:),'r.');
    
    % For every pixel
    % for every pixel want to figure out which lattice site it is
    for ii=1:numel(z3)   
        m1 = M(1,ii)-n1i+1;
        m2 = M(2,ii)-n2i+1;        
        ind=(m1-1)*numel(n2)+m2;              
        z3(ii) = ind;
        z3(ii) = mod(ind,2);
    end
    
    subplot(122);
    imagesc(x2,y2,z3);
    hold on
    plot(x_site(:),y_site(:),'r.');
    
end

%%

[nn1,nn2]=meshgrid(n1,n2);

N=[nn1(:)' ; nn2(:)'];                      % All points
P = repmat(p,[1 length(N)]);        % Phase vector        


Rn=A*(N+P);                                 % Positino of every lattice site

Xn = reshape(Rn(1,:),size(Zbin));
Yn = reshape(Rn(2,:),size(Zbin));

out = struct;
out.p = [opts.p1 opts.p2];
out.a1 = opts.a1;
out.a2 = opts.a2;
out.n1 = n1;
out.n2 = n2;
out.Zbin = Zbin;

p1 = opts.p1;
p2 = opts.p2;
a1 = opts.a1;
a2 = opts.a2;
out.site2px = @(n1,n2) (n1+p1)*a1 + (n2+p2)*a2;

% Lattice Spacing in pixels
lattice_spacing_px = mean([norm(a1) ...
    norm(a2)]);

%Lattice spacing in um
lattice_spacing_um = 0.527; 

out.lattice_spacing_um = lattice_spacing_um;
out.lattice_spacing_px = lattice_spacing_px;
out.site2um = @(n1,n2) ((n1+p1)*a1 + (n2+p2)*a2)*(lattice_spacing_um/lattice_spacing_px);
%% Threshold Detection


Zest=Zbin(:);
Zest(Zest==0)=[];
cluster_number=2;
[idx,c,sumD,D]=kmeans(Zest,cluster_number);

[c,inds]=sort(c,'ascend');
sumD=sumD(inds);

dd = zeros(numel(c),1);
thresh = zeros(numel(c),1);
for kk=1:numel(c)
    ind = inds(kk);
    nThis = sum(idx==ind);
    dd(kk) = sqrt(sumD(kk)/nThis);
    thresh(kk) = round(max(Zest(idx==ind)));
end
thresh(end)=[];

cluster_overlap = zeros(numel(thresh),1);

        z2 = data(data>thresh(1));
        % The n=1 (atom) distribution is fit with a simple gaussian
        % distribution
        % pd1 = fitdist(z2,'normal'); % This works too, but use MLE to make it
        % the same output
        [pdf1_c,pdf1_cint] = mle(z2,'distribution','normal');
        
        I0=pdf1_c(1); % center and then the 
        



for kk=1:numel(thresh)
    cluster_overlap(kk)=dd(kk);
end

out.ClusterNumber = cluster_number;
out.ClusterThreshold = thresh;
out.ClusterCentroid = c;
out.ClusterRadius = dd;

end
