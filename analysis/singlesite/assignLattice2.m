function out=assignLattice2(x,y,z,opts)
% Author : C. Fujiwara
%
% Given a fluoresence image input (z) with corresponding pixel position
% cooridinates (x,y), this function assigns fluoresnce to each lattice
% site.
%
% The phase is found my optizing the total contrast of the image which is
% gotten from the 

tic;
fprintf('lattice binning ...');

if nargin < 4
    opts = struct;
    opts.a1 = [1.3316; 2.3134];
    opts.a2 = [-2.2851; 1.3896];
    opts.ScaleFactor = 5;
end

if ~isfield(opts,'ScaleFactor')
    opts.ScaleFactor = 5;
end

if ~isfield(opts,'phase')
    opts.FindPhase = 1;
    Ncoarse = 4;            % Ncoarse x Ncoasre search phase
    Nfine = 20;             % 2 x Nfine search for phase
else 
    opts.FindPhase = 0;
end

% Matrix which defines lattice basis
A = [opts.a1 opts.a2];                  % lattice vectors

%% Image preparation

% Rescale the image to avoid rounding errors when moving pixels
z2 = imresize(z,opts.ScaleFactor)/(opts.ScaleFactor^2); % scale amplitude to preserve norm
x2 = linspace(x(1),x(end),size(z2,2));
y2 = linspace(y(1),y(end),size(z2,1));

% Create 2xN vectors of all pixels positions
[X,Y]=meshgrid(x2,y2);                  % matrix of X and Y
R=[X(:)' ; Y(:)'];                      % All points

% 1xN vector all counts per pixel
Z = z2(:);                              % counts per re-scaled pixel

% Solve for all position in terms of the lattice basis
N0 = inv(A)*R;

%% Find bounds on lattice indeces
% Given the lattice basis and the image bounds, determine some reasonable
% bounds on the extremal lattice indeces.

% Solve for positionof each corner in terms of the lattice basis
n1 = inv(A)*[x(1); y(1)];
n2 = inv(A)*[x(1); y(end)];
n3 = inv(A)*[x(end); y(1)];
n4 = inv(A)*[x(end); y(end)];

% Collect the four corners
Nc = [n1 n2 n3 n4];

% Bounds and vector for lattice basis vector 1
i1 = floor(min(Nc(1,:)));
i2 = ceil(max(Nc(1,:)));
ivec = i1:i2;

% Bounds and vector for lattice basis vector 2
j1 = floor(min(Nc(2,:)));
j2 = ceil(max(Nc(2,:)));
jvec = j1:j2;

% Matrix of all lattice positions
[imat,jmat] = meshgrid(ivec,jvec);


%% Calculate the Phase

% CF isn't sure why this works, BUT WHATEVETA FOURIER TRANSFORMS
p1=mod(1-(angle(sum(exp(-1i*2*pi*(opts.k1(1).*X+opts.k1(2).*Y)).*z2,'all')))/(2*pi),1);
p2=mod(1-(angle(sum(exp(-1i*2*pi*(opts.k2(1).*X+opts.k2(2).*Y)).*z2,'all')))/(2*pi),1);

%% Final Binning
fprintf('binning ...');

p = [p1; p2];                       % Phase        
P = repmat(p,[1 length(Z)]);        % Phase vector        
M = round(N0-P);                    % round pixel to lattice site
ibad = logical((M(1,:)<i1) + (M(1,:)>i2) + (M(2,:)>j2) + (M(2,:)<j1));

M(:,ibad) = [];

Zbin = zeros(j2-j1+1,i2-i1+1);
for ii=1:size(M,2)    
    m1 = M(1,ii)-i1+1;
    m2 = M(2,ii)-j1+1;
    Zbin(m2,m1) = Zbin(m2,m1) + Z(ii);
end

Zall = Zbin;
Zall = Zall(:);
disp([p1 p2]);
%%
[nn1,nn2]=meshgrid(ivec,jvec);

N=[nn1(:)' ; nn2(:)'];                      % All points
P = repmat(p,[1 length(N)]);        % Phase vector        


Rn=A*(N+P);                                 % Positino of every lattice site

out = struct;
out.p = [p1 p2];
out.a1 = opts.a1;
out.a2 = opts.a2;
out.n1 = ivec;
out.n2 = jvec;
out.Zbin = Zbin;
out.R = Rn;

t2 = toc;
disp(['done (' num2str(t2) ' seconds)']);


end