function out = findLatticeK(fx,fy,Zf,opts)

    pixelsize = 16/81; % um/px
    a = .527;          % Lattice spacing in um
    a = a/pixelsize;   % Lattice spacing in px
    
    if nargin ==4
       if isfield(opts,'doScale') && isfield(opts,'ScaleFactor') && opts.ScaleFactor > 1 && opts.doScale
          a = a*opts.ScaleFactor; 
       end
    end
    


% Lattice wavevector
kL = 1/a;

[fxx,fyy]=meshgrid(fx,fy);
fmat = sqrt(fxx.^2+fyy.^2);

% Make Mask about the lattice k-vectors
fmat(fmat>(kL*1.1)) = 0;
fmat(fmat<(kL*0.9)) = 0;

% Mask the data
ss = (fmat~=0).*abs(Zf);

% Find maximum peak in quadrant 1
tt=atan2(fyy,fxx);
ii = (tt<(0.25*pi)).*(tt>(-0.25*pi));

A=(fxx>0).*(fyy>0).*ss;

A=ii.*ss;

[~,I] = max(A,[],"all","linear");
[dim1, dim2] = ind2sub(size(A),I);
k1 = [fx(dim2); fy(dim1)];

% Pixel radius to fit about
pxR = 10;

% 2D Gaussian Function (to fit peak)
myfit = fittype('A*exp(-(xx-xc).^2/(2*s^2)).*exp(-(yy-yc).^2/(2*s^2))+b',...
    'independent',{'xx','yy'},'coefficients',{'A','xc','yc','s','b'});
fitopt = fitoptions(myfit);
    
% Fit the peak in the Quadrant One
tic
fxs = fx(dim2+[-pxR:pxR]);
fys = fy(dim1+[-pxR:pxR]);    
[fxxx,fyyy]=meshgrid(fxs,fys);
zsub = Zf(dim1+[-pxR:pxR],dim2+[-pxR:pxR]);
zz = abs(zsub);

ampl = max(max(zz)) - min(min(zz));
bkgd =min(min(zz));

fitopt.StartPoint = [ampl k1(1) k1(2) .001 bkgd];    

x = fxxx(:);
y = fyyy(:);
z = zz(:);

ilow = [(z-bkgd)<(ampl/2)];
x(ilow) = [];
y(ilow) = [];
z(ilow) = []; 

[fout1,gof1,output1] = fit([fxxx(:),fyyy(:)],zz(:),myfit,fitopt);    
% [fout1a,gof1a,output1a] = fit([x,y],z,myfit,fitopt);    



k1 = [fout1.xc; fout1.yc];
s1 = fout1.s;
t2=toc;
        
% Find peak in quadrant 2
A=(fxx>0).*(fyy<0).*ss;
    ii = (tt<(0.75*pi)).*(tt>(0.25*pi));

A = ii.*ss;

[~,I] = max(A,[],"all","linear");
[dim1, dim2] = ind2sub(size(A),I);
k2 = [fx(dim2); fx(dim1)];
    
% Fit the peak in the Quadrant One
tic
fxs = fx(dim2+[-pxR:pxR]);
fys = fx(dim1+[-pxR:pxR]);    
[fxxx,fyyy]=meshgrid(fxs,fys);
zsub = Zf(dim1+[-pxR:pxR],dim2+[-pxR:pxR]);
zz = abs(zsub);
fitopt.StartPoint = [max(max(zz)) k2(1) k2(2) .001 min(min(abs(zsub)))];    
fout2 = fit([fxxx(:),fyyy(:)],zz(:),myfit,fitopt);    
k2 = [fout2.xc; fout2.yc];   
s2 = fout2.s;
t2=toc;

%% Output

out = struct;
out.k1 = k1;
out.k2 = k2;
out.s1 = s1;
out.s2 = s2;
out.Fit1 = fout1;
out.Fit2 = fout2;
 
end

