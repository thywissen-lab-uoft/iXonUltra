function out = findLatticeK(fx,fy,Zf)

pixelsize = 16/81; % um/px

a = .527;          % Lattice spacing in um
a = a/pixelsize;   % Lattice spacing in px

% Lattice wavevector
kL = 1/a;

[fxx,fyy]=meshgrid(fx,fy);
fmat = sqrt(fxx.^2+fyy.^2);

% Make Mask about the lattice k-vectors
fmat(fmat>(kL*1.05)) = 0;
fmat(fmat<(kL*0.95)) = 0;

% Mask the data
ss = (fmat~=0).*abs(Zf);

% Find maximum peak in quadrant 1
A=(fxx>0).*(fyy>0).*ss;
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
fitopt.StartPoint = [max(max(zz)) k1(1) k1(2) .01 min(min(abs(zsub)))];    
fout1 = fit([fxxx(:),fyyy(:)],zz(:),myfit,fitopt);    
k1 = [fout1.xc; fout1.yc];   
s1 = fout1.s;
t2=toc;
        
% Find peak in quadrant 2
A=(fxx>0).*(fyy<0).*ss;
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
fitopt.StartPoint = [max(max(zz)) k2(1) k2(2) .01 min(min(abs(zsub)))];    
fout2 = fit([fxxx(:),fyyy(:)],zz(:),myfit,fitopt);    
k2 = [fout2.xc; fout2.yc];   
s2 = fout2.s;
t2=toc;

out = struct;
out.k1 = k1;
out.k2 = k2;
out.s1 = s1;
out.s2 = s2;
out.Fit1 = fout1;
out.Fit2 = fout2;
 
% disp(['k1 ' num2str(k1(1)) ',' num2str(k1(2))]);
% disp(['k2 ' num2str(k2(1)) ',' num2str(k2(2))]);

end

