function output = ixon_getBoxData(ixondata,xVar)
%GETERFDATA Summary of this function goes here
%   Detailed explanation goes here

if nargin == 1
   xVar = 'ExecutionDate'; 
end

%% Sort the data by the parameter given
params=[ixondata.Params];
X=[params.(xVar)];

[X,inds]=sort(X,'ascend');
ixondata=ixondata(inds);

% Make sure its Nx1
X = reshape(X,[length(X) 1]);

%% Grab the Erf Fit outputs

bgROI=ixondata(1).BoxCount.bgROI;        % ROI for calculating bgkd

for kk=1:length(ixondata)
   for nn=1:length(ixondata(kk).BoxCount)
        bc=ixondata(kk).BoxCount(nn);               % Grab the fit   
        N(kk,nn)=bc.Ncounts;        % Number of counts (w/ bkgd removed)
        Nraw(kk,nn)=bc.Nraw;          % Raw of number of counts
        Nbg(kk,nn)=bc.Nbkgd;          % Bakcground number of counts
        nbg(kk,nn)=bc.nbkgd;          % Background counts/px
        Xc(kk,nn)=bc.Xc;              % X center of mass
        Yc(kk,nn)=bc.Yc;              % Y center of mass
        Xs(kk,nn)=bc.Xs;              % X standard deviation
        Ys(kk,nn)=bc.Ys;              % Y standard deviation                      
        Zs(kk,nn)=bc.Xs;              % Y standard deviation      
   end        
end

output = struct;
output.FileNames    = {ixondata.Name}';
output.xVar         = xVar;
output.X            = X;
output.Params       = params;
output.Units        = [ixondata.Units];
output.Flags        = [ixondata.Flags];
output.FitType      = 'box';
output.PixelSize    = ixondata(1).PixelSize;
output.Magnification = ixondata(1).Magnification;

% Assign fit outputs
output.bgROI        = bgROI;
output.N            = N;
output.Xc           = Xc;
output.Yc           = Yc;
output.Xs           = Xs;
output.Ys           = Ys;
output.Zs           = Zs;
output.nbg          = nbg;

end

