function output = ixon_getGaussData(ixondata,xVar)
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

for kk=1:length(ixondata)
    for ll = 1:size(ixondata(kk).GaussFit,1)
       for nn=1:length(ixondata(kk).GaussFit(ll,:))            
            fout=ixondata(kk).GaussFit{ll,nn};             % Grab the fit             
            Xc(kk,ll,nn)=fout.Xc;Yc(kk,ll,nn)=fout.Yc;        % X and Y center
            Xs(kk,ll,nn)=fout.Xs;Ys(kk,ll,nn)=fout.Ys;        % X and Y sigma
            A(kk,ll,nn)=fout.A;                            % Amplitude
            nbg(kk,ll,nn)=fout.nbg;                        % Background
            N(kk,ll,nn)=2*pi*Xs(kk,ll,nn)*Ys(kk,ll,nn)*A(kk,ll,nn);  % Number of counts
       end 
   end        
end

output = struct;
output.FileNames    = {ixondata.Name}';
output.xVar         = xVar;
output.X            = X;
output.Params       = params;
output.Units        = [ixondata.Units];
output.Flags        = [ixondata.Flags];
output.FitType      = 'gauss';
output.PixelSize    = ixondata(1).PixelSize;
output.Magnification = ixondata(1).Magnification;

% Assign fit outputs
output.N            = N;
output.Xc           = Xc;
output.Yc           = Yc;
output.Xs           = Xs;
output.Ys           = Ys;
output.nbg          = nbg;


end

