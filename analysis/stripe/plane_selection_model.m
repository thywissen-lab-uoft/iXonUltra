function out = plane_selection_model(opts)
out = struct;
if nargin ==1
    opts = struct;
end
% Feshbach field in Gauss; assumed in the z direction
Bfb = 120;

% Calculate field gradient using selection frequency
G_freq = 45;            % vertical QP gradient in kHz/plane
G_freq = 90;            % vertical QP gradient in kHz/plane
zeman = 2600;           % 2500 kHz/Gauss for transition
a = 532e-9*1e2;         % plane separation in cm
G = (G_freq/zeman)/a;   % Field gradient in Gauss/cm
% G = G*1.05;

% Vertical QP gradient G/cm
% G = 210;

% Vertical magnetic field per plane from the QP gradient
dB = G*a;

% QP coils center (in um)
% atomic center is assumed to be at x=z=0;
x0 = 240;0;
z0 = 40;

% Shim Field (in gauss)
Bshim_x =-3.9;-7.7;0;
Bshim_x = -x0*1e-4*G/2;

Bshim_z = 0;

out.Bfb = Bfb;
out.G = G;
out.dB = dB;
out.x0 = x0;
out.z0 = z0;
out.Bshim_x = Bshim_x;
out.Bshim_z = Bshim_z;
%% String Stuff

% Critical radius at which you cross into another plane
rC =  r_critical(Bfb,G);

% Summary Text
str = ['$B_\mathrm{fb0}:' num2str(Bfb) '~\mathrm{G}$' newline  ... 
    '$G_0:' num2str(round(G)) '~\mathrm{G/cm}, (' num2str(round(dB*1e3,1)) ...
    '~\mathrm{mG},' num2str(round(dB*zeman,1)) '~\mathrm{kHz})/\mathrm{plane}$' newline ...
    '$(x_0,z_0) : (' num2str(x0) ',' num2str(z0) ')~\mu\mathrm{m}$' newline ...
    '$r_c:' num2str(round(rC*1e4)) '~\mu\mathrm{m}$' newline ...
    '$B_\mathrm{s0}: (' num2str(Bshim_x) ',' num2str(Bshim_z) ')~\mathrm{G}$'];
%% Field Functions

% Feshbach field function is assumed to be ideal constant field IDEAL
    function [Br,Bz]=B_feshbach(B0,r,z)% positions in cm
        Br = zeros(size(r,1),size(r,2));
        Bz = B0*ones(size(r,1),size(r,2));
    end

% Quadrulpolar field % positions in cm IDEAL
    function [Bx,Bz]=B_qp(G,x0,z0,x,z) 
        Bx = -G/2.*(x-x0);
        Bz = G.*(z-z0);
    end

% Shim field field % positions in cm IDEAL
    function [Bx,Bz]=B_shim(Bsx,Bsz,x,z)
        Bx = Bsx;
        Bz = Bsz;
    end

    function rC = r_critical(B0,G)
        rC=sqrt(2*G*a/B0/(G/(2*B0))^2);
    end

%% Calculate Total Field

% Horizontal and Vertical Vector [um]
r = linspace(-200,200,1e3);
x = linspace(-2,2,1e3);
[R,Z]=meshgrid(r,x);

% Calculate individual fields
[Bx_qp,Bz_qp] = B_qp(G,x0*1e-4,z0*1e-4,1e-4*R,1e-4*Z);
[Bx_fb,Bz_fb] = B_feshbach(Bfb,1e-4*R,1e-4*Z);
[Bx_s,Bz_s] = B_shim(Bshim_x,Bshim_z);

Bx = Bx_qp + Bx_fb+Bx_s;
Bz = Bz_qp + Bz_fb+Bz_s;
B = (Bx.^2+Bz.^2).^(1/2);



[Bx_qp0,Bz_qp0] = B_qp(G,x0*1e-4,z0*1e-4,0,0);
[Bx_fb0,Bz_fb0] = B_feshbach(Bfb,0,0);
[Bx_s0,Bz_s0] = B_shim(Bshim_x,Bshim_z);
Bx0 = Bx_qp0 + Bx_fb0 + Bx_s0;
Bz0 = Bz_qp0 + Bz_fb0+ Bz_s0;
B0 = (Bx0.^2+Bz0.^2).^(1/2);

out.B0 = B0;
%% Plotting

hf=figure(10);
hf.Color='w';
hf.Position=[100 100 1200 500];
clf

% Contour Plot
subplot(1,4,[1 2]);
Bc = Bfb;
Bc = 118.7;
Bc = B0;
levels = Bc +[-15:1:15]*dB;
[c,h]=contour(R,Z,B,'levels',dB);
h.LevelList = levels;
h.LineWidth = 2;
h.ShowText = ['on'];
% title(str,'interpreter','latex');
text(.01,.01,str,'units','normalized','interpreter','latex',...
    'verticalalignment','bottom','horizontalalignment','left','fontsize',10,...
    'backgroundcolor',[1 1 1 .85])
xlabel('horizontal position (\mum)');
ylabel('axial position (\mum)');
set(gca,'fontsize',10,'fontname','times','box','on',...
    'linewidth',1)
grid on
title('magnetic field profile');

%% Horizontal Cut Plut
r= linspace(-200,200,1e4);

[Bx_qp_cut,Bz_qp_cut] = B_qp(G,x0*1e-4,z0*1e-4,1e-4*r,0);
[Bx_fb_cut,Bz_fb_cut] = B_feshbach(Bfb,1e-4*r,0);
[Bx_s_cut,Bz_s_cut] = B_shim(Bshim_x,Bshim_z);
Bx_cut = Bx_qp_cut + Bx_fb_cut + Bx_s_cut;
Bz_cut = Bz_qp_cut + Bz_fb_cut+ Bz_s_cut;
B_cut = (Bx_cut.^2+Bz_cut.^2).^(1/2);


dB_relative = (B_cut-B0)/dB;
[val,ind]=min(abs(dB_relative));
i2 = find(abs(dB_relative(ind:end))>1,1);
rC = r(i2+ind-1);
out.rC = rC;

subplot(2,4,3);
plot(r,dB_relative,'linewidth',2);
hold on
xlabel('horizontal position (\mum)');
ylabel('$\Delta B/\Delta B_\mathrm{plane}$ (G)','interpreter','latex')
set(gca,'fontsize',10,'fontname','times','box','on',...
    'linewidth',1)
grid on
title('horizontal field profile (z=0)');
pR=plot(rC,1,'ro','linewidth',1,'markerfacecolor','r');
legend(pR,[num2str(round(rC)) 'um'])

%% Field Dependece at x=0,z=0

Bfb_vec = linspace(.95,1.05,1e4)*Bfb;
G_vec = linspace(.1,2,1e4)*G;
Bsx_vec = Bshim_x+linspace(-5,5,1e4);

% Variable fb current
[Bx_qp_chfb,Bz_qp_chfb] = B_qp(G,x0*1e-4,z0*1e-4,0,0);
[Bx_fb_chfb,Bz_fb_chfb] = B_feshbach(Bfb_vec,0,0);
[Bx_s_chfb,Bz_s_chfb] = B_shim(Bshim_x,Bshim_z);
Bx_chfb = Bx_qp_chfb + Bx_fb_chfb + Bx_s_chfb;
Bz_chfb = Bz_qp_chfb + Bz_fb_chfb+ Bz_s_chfb;
B_chfb = (Bx_chfb.^2+Bz_chfb.^2).^(1/2);

% Variable Gradient
[Bx_qp_chG,Bz_qp_chG] = B_qp(G_vec,x0*1e-4,z0*1e-4,0,0);
[Bx_fb_chG,Bz_fb_chG] = B_feshbach(Bfb,0,0);
[Bx_s_chG,Bz_s_chG] = B_shim(Bshim_x,Bshim_z);
Bx_chG = Bx_qp_chG + Bx_fb_chG + Bx_s_chG;
Bz_chG = Bz_qp_chG + Bz_fb_chG+ Bz_s_chG;
B_chG = (Bx_chG.^2+Bz_chG.^2).^(1/2);

% Variable Shim
[Bx_qp_chS,Bz_qp_chS] = B_qp(G,x0*1e-4,z0*1e-4,0,0);
[Bx_fb_chS,Bz_fb_chS] = B_feshbach(Bfb,0,0);
[Bx_s_chS,Bz_s_chS] = B_shim(Bsx_vec,Bshim_z);
Bx_chS = Bx_qp_chS + Bx_fb_chS + Bx_s_chS;
Bz_chS = Bz_qp_chS + Bz_fb_chS+ Bz_s_chS;
B_chS = (Bx_chS.^2+Bz_chS.^2).^(1/2);

% Shim Sensitivity
subplot(2,4,4);
dB_relative = (B_chS-B0)/dB;
[~,ind]=min(abs(dB_relative));
i2 = find(abs(dB_relative(ind:end))>1,1);
Bsx_c = Bsx_vec(i2+ind-1);
plot(Bsx_vec,(B_chS-B0)/dB,'linewidth',2);
xlim([min(Bsx_vec) max(Bsx_vec)]);
xlabel('shim $B_{sx}$ (G)','interpreter','latex')
ylabel('$\Delta B/\Delta B_\mathrm{plane}$','interpreter','latex')
ylim(0+[-1.5 1.5]);
grid on
set(gca,'fontsize',10,'box','on','linewidth',1,'fontname','times')
title('shim sensitivity (r=z=0)');
hold on
pR=plot(Bsx_c,dB_relative(i2+ind-1),'ro','linewidth',1,'markerfacecolor','r');
dB_shim  = Bsx_c-Bshim_x;
legend(pR,['\Delta B : ' num2str(round(dB_shim,2)) ' G/plane'],'location','best')
out.dShim_G = dB_shim;

% Gradient Sensitivity
subplot(2,4,7);
dB_relative = (B_chG-B0)/dB;
[~,ind]=min(abs(dB_relative));
i2 = find(abs(dB_relative(ind:end))>1,1);
G_c = G_vec(i2+ind-1);
plot(G_vec/G,(B_chG-B0)/dB,'linewidth',2);
xlabel('gradient $G/G_0$','interpreter','latex')
ylabel('$\Delta B/\Delta B_\mathrm{plane}$','interpreter','latex')
ylim(0+[-1.5 1.5]);
grid on
set(gca,'fontsize',10,'box','on','linewidth',1,'fontname','times')
title('gradient sensitivity (r=z=0)');
hold on
pR=plot(G_c/G,dB_relative(i2+ind-1),'ro','linewidth',1,'markerfacecolor','r');
ppm = (G_c/G-1)*1e6;
legend(pR,[num2str(round(ppm)) ' ppm/plane'],'location','best')
out.dG_ppm = ppm;

% Feshbach Sensitivity
subplot(2,4,8);
dB_relative = (B_chfb-B0)/dB;
[~,ind]=min(abs(dB_relative));
i2 = find(abs(dB_relative(ind:end))>1,1);
Bfb_c = Bfb_vec(i2+ind-1);
plot((Bfb_vec)/Bfb,(B_chfb-B0)/dB,'linewidth',2);
ylim(0+[-1.5 1.5]);
xlabel('feshbach $B_\mathrm{fb}/B_\mathrm{fb0}$','interpreter','latex')
ylabel('$\Delta B/\Delta B_\mathrm{plane}$','interpreter','latex')
grid on
set(gca,'fontsize',10,'box','on','linewidth',1,'fontname','times')
title('feshbach sensitivity (r=z=0)');
hold on
pR=plot(Bfb_c/Bfb,dB_relative(i2+ind-1),'ro','linewidth',1,'markerfacecolor','r');

ppm = (Bfb_c/Bfb-1)*1e6;
out.dFB_ppm = ppm;

legend(pR,[num2str(round(ppm)) ' ppm/plane'],'location','best')

end

