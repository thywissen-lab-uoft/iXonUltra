
% useful constants
m = 39.96399848*1.66053906660*10^-27;
w_xdt_x = 2*pi*66.8;
w_xdt_y = 2*pi*59.7; % old
h = 6.62607*10^-34;
hbar = h/(2*pi);
kB = 1.380649*10^-23;

% size of the cloud to temp
xs_um = [digdata.Xs_um];
ys_um = [digdata.Ys_um];

T_x = m*w_xdt_x^2*(xs_um*1e-6).^2 / kB *1e9;
T_y = m*w_xdt_y^2*(ys_um*1e-6).^2 / kB *1e9;


% peak filling
latt_depth = [digdata.X];
X = latt_depth;
Ux = unique(latt_depth);

Z_0 = 0;
fill_err = 0;
for kk = 1:length(Ux)
   if opts.ForceAverage         
        inds = 1:length(X);
   else
        x = Ux(kk);
        inds = find([X==x]);
   end

   nImg = length(inds);                          % Number of images
   Z = mean(digdata.Zdig(:,:,inds),3);          % Average Image for this X variable
   Zr = mean(digdata.Zr(:,inds),2);             % Average Radial data for this X variable
   err = std(digdata.Zr(:,inds),1,2)/sqrt(nImg);  % Standard error at each radius 
   r = mean(digdata.r(:,inds),2);               % Radial vector
   
   inds = find(0.5>Zr & Zr>0);
   Zr = Zr(inds);
   Z_0(end+1) = mean(Zr(1:3));
   fill_err(end+1) = std(Zr(1:3));
   
   disp(Zr)
   
end

Z_0(1) = [];
fill_err(1) = [];
npeak = Z_0;

% fit first few data points to linear
p = polyfit(latt_depth(1:end-2),T_x(1:end-2),1);



% Plot against lattice depth
co=get(gca,'colororder');

hf = figure('color','w');


subplot(2,1,1)
plot(latt_depth,T_x,'o','color',co(1,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(1,:),'markeredgecolor',co(1,:)*.5);
hold on
xx = [-0.5:0.1:10];
plot(xx,polyval(p,xx));
xlabel(digdata.xVar,'Interpreter','none')
ylabel('Gauss. X temp. (nK)')

legend('Data',['y ='+string(p(1))+'x + '+string(p(2))])

subplot(2,1,2)
errorbar(Ux,npeak,fill_err,'o','color',co(1,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(1,:),'markeredgecolor',co(1,:)*.5);
xlabel(digdata.xVar,'Interpreter','none')
ylabel('Peak charge filling')






