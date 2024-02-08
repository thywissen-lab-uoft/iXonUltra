function out=binLattice(x,y,z,opts)
% Author : C. Fujiwara
%
% Given a fluoresence image input (z) with corresponding pixel position
% cooridinates (x,y), this function assigns fluoresnce to each lattice
% site.
%
% The phase is found my optizing the total contrast of the image which is
% gotten from the 

% Matrix which defines lattice basis
A = [opts.a1 opts.a2];                  % lattice vectors

%% Image preparation

% Rescale the image to avoid rounding errors when moving pixels
z2 = imresize(z,opts.ScaleFactor,'method','bilinear')/(opts.ScaleFactor^2); % scale amplitude to preserve norm
x2 = linspace(x(1),x(end),size(z2,2));
y2 = linspace(y(1),y(end),size(z2,1));


% Create 2xN vectors of all pixels positions
[X,Y]=meshgrid(x2,y2);                  % matrix of X and Y
R=[X(:)' ; Y(:)'];                      % All points

if isfield(opts,'PixelThreshold') && isfield(opts,'UsePixelThreshold')
    z2(z2<(opts.PixelThreshold/opts.ScaleFactor^2))=0;    
end

% 1xN vector all counts per pixel
Z = z2(:);                              % counts per re-scaled pixel

% Solve for all position in terms of the lattice basis
N0 = inv(A)*R;

p = [opts.p1; opts.p2];                       % Phase        
P = repmat(p,[1 length(Z)]);        % Phase vector        
M = round(N0-P);                    % round pixel to lattice site

% aa = reshape(M(1,:),[size(X,1) size(X,2)]);
% bb = reshape(M(2,:),[size(X,1) size(X,2)]);
% cc = mod(bb+aa,2);
% 
% dd = z2.*cc;
% keyboard
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

% Remove points outside of the desired ROI.
ibad = logical((M(1,:)<n1i) + (M(1,:)>n1f) + (M(2,:)>n2f) + (M(2,:)<n2i));
M(:,ibad) = [];
Z(ibad) = [];

% CF I forget what this part of the code does with M0
% M0 = N0-P;
% M0(:,ibad)=[];

Zbin = zeros(n2f-n2i+1,n1f-n1i+1);
Znum = zeros(n2f-n2i+1,n1f-n1i+1);
for ii=1:size(M,2)    
    m1 = M(1,ii)-n1i+1;
    m2 = M(2,ii)-n2i+1;
    
    Znum(m2,m1) = Znum(m2,m1)+1;
    Zbin(m2,m1) = Zbin(m2,m1) + Z(ii);  
end

Zall = Zbin;
Zall = Zall(:);
Zbin2 = Zbin;
Zpercent = Znum/max(Znum,[],'all');
Zbin2(Zpercent<.9)=NaN;

Zbin = Zbin2;


% [cc,rr]=meshgrid(1:size(Zbin,1),size(Zbin,2))
% min(cc(isnan(Zbin)))
% m




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

out.R0 = [Xn(1) Yn(1)];     % Position of first lattice site in pixels
out.N0 = [N(1,1) N(2,1)];   % Value of n1 n2 corresponding to the first lattice site

% Dont store these for space considerations
% out.Xn = Xn;
% out.Yn = Yn;


% 
% if isfield(opts, 'DigitizationThreshold')
%     
%    out.Zdig = Zbin>=opts.DigitizationThreshold; 
%    
%   
% 
% end

%% Digital Box Count

%   
%             Natoms = sum(sum(Zdig));        % Total number of atoms
%             Nsites = length(Zdig(:));       % Total number of sites        
% %             filling = Natoms/Nsites;        % Average filling fraction
%             
%   
% 
%             % Calculate center of mass
%             Xc=sum(zX.*x);
%             Yc=sum(zY.*y);          
% 
%             % Calculate central second moment/variance and the standard
%             % deviation
%             X2=sum(zX.*(x-Xc).^2); % x variance
%             Xs=sqrt(X2); % standard deviation X
%             Y2=sum(zY.*(y-Yc).^2); % x variance
%             Ys=sqrt(Y2); % standard deviation Y               
% 
%             LatticeDig(k).Natoms = Natoms;
% %             LatticeDig(k).Filling = filling;
%             LatticeDig(k).Xc = Xc;
%             LatticeDig(k).Yc = Yc;
%             LatticeDig(k).Xs = Xs;
%             LatticeDig(k).Ys = Ys;   

end
