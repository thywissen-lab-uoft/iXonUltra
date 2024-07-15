function out=binLattice2(x,y,z,opts)
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

p = [opts.p1; opts.p2];             % Phase        

% Image pixel size in units of camera pixels
% dX = x(2)-x(1);
% dY = y(2)-y(1);


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

%% Image preparation

% Rescale the image to avoid rounding errors when moving pixels
% scale amplitude to preserve norm
z2 = imresize(z,opts.ScaleFactor,'method','bilinear')/(opts.ScaleFactor^2); 
x2 = linspace(x(1),x(end),size(z2,2));
y2 = linspace(y(1),y(end),size(z2,1));

[xx,yy]=meshgrid(x2,y2);
%%  Iterate over every pixel

Zbin = zeros(n2f-n2i+1,n1f-n1i+1);
Znum = zeros(n2f-n2i+1,n1f-n1i+1);
disp('iterarting');
for kk=1:numel(z2)
%     tic
    r = [xx(kk);yy(kk)];
    n0 = inv(A)*r;
    n  = round(n0-p);

    row = n(1)-n1i+1;
    col = n(2)-n2i+1;
    
    
%     try
        Znum(row,col) = Znum(row,col)+1;
        Zbin(row,col) = Zbin(row,col) + z2(kk);  
%     catch ME
%         warning('uhoh');
%     end
%     toc
    
end

disp('done');
%%

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

%%

end
