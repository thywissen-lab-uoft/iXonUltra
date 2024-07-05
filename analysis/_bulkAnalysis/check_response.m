%% Constants
h = 6.626e-34;
hbar = h/(2*pi);
kb = 1.380649e-23;
aL = 527e-9;
amu = 1.660538921e-27;
m = 39.964008*amu;

%% Define the lookup tables for R and E
global Rvalues;
Rvalues_unscaled = table2array(readtable('Rvalues_unscaled_64_4Hz_200.csv'));
Rvalues = -1j*(aL/pi)*Rvalues_unscaled;

global energies;
energies_Hz = importdata('EnergyHz_64_4Hz_200.txt');
energies = h*(energies_Hz);

%% Trap parameters
wXDT = 2*pi*42.5; %2*pi*Hz
T = 80e-9; %K
G = 2*pi*35; %2*pi*Hz
amp_desired = 0.8; %um

f = [20:1:150];

%% Create data points

d = [];
for ff = 1:length(f)
    d(ff) = hbar*2*pi*f(ff)*amp_desired/(aL^2*m*wXDT^2*sqrt(qfit_real(T,G,2*pi*f(ff))^2+qfit_imag(T,G,2*pi*f(ff))^2))/3.55; %V
end

%% Plot Anticipated response
f33 = figure;
f33.Color='w';
clf(f33);
plot(f,d,'k.-')
hold on;
xlabel('drive frequency (Hz)');
ylabel('drive amplitude (V)');

% plot(xx,yy,'r')
% plot(x2,y2,'b')
xlim([20,130])




%% Fit the response


% plot(f,polyval(pp,f-f0),'r-');


[d0,ind]=min(d);
f0 = f(ind);
pp = polyfit(f-f0,d,4);

WH = [f>=f0]';
WL = [f<=f0]';

myfit = fittype(@(a8,a6,a4,a2,x) d0 + a2*(x-f0).^2 + a4*(x-f0).^4+ a6*(x-f0).^6+ a8*(x-f0).^8,...
    'independent',{'x'},'coefficients',{'a8','a6','a4','a2'});
opt = fitoptions(myfit);
opt.StartPoint = [0 0 pp(1) pp(3)];
opt.Weights = WH;
opt.Robust = 'bisquare';

foutH = fit(f',d',myfit,opt);

opt.Weights = WL;
foutL = fit(f',d',myfit,opt);

plot(f(f>=f0),feval(foutH,f(f>=f0)),'r-','linewidth',1)
plot(f(f<=f0),feval(foutL,f(f<=f0)),'b-','linewidth',1)

s = ['$y = y_0 + a_2(x-x_0)^2 + a_4(x-x_0)^4 + a_6(x-x_0)^6 + a_8(x-x_0)^8$'];
s = [s newline '$(x_0,y_0) = ' num2str(round(f0,3)) ',' num2str(round(d0,4)) '$'];
s = [s newline '$(a_2,a_4,a_6,a_8)_L = ($' num2str(foutL.a2,'%.2e') ',' ...
    num2str(foutL.a4,'%.2e') ', ' ...
    num2str(foutL.a6,'%.2e') ', ' ...
    num2str(foutL.a8,'%.2e') '$)$'];
s = [s newline '$(a_2,a_4,a_6,a_8)_H = ($' num2str(foutH.a2,'%.2e') ',' ...
    num2str(foutH.a4,'%.2e') ', ' ...
    num2str(foutH.a6,'%.2e') ', ' ...
    num2str(foutH.a8,'%.2e') '$)$'];

text(.01,.98,s,'interpreter','latex','units','normalized','fontsize',12,...
    'verticalalignment','top');


% plot([20 40], [1.63,1.63],'r')
% plot([98 100], [4,4],'r')
% xlabel('Drive Frequency (Hz)')
% ylabel('Drive Amplitude (V)')
% text(22,3.4,'y =a(x-b)^2 + c','color','r')
% text(22,3.2,['a = ' num2str(fout.a)],'color','r')
% text(22,3,['b = ' num2str(fout.b)],'color','r')
% text(22,2.8,['c = ' num2str(fout.c)],'color','r')
% text(22,2.6,'y= -0.0084x + 0.9340','color','b')

% %% Fit the response function 
% myfunc = @(A,x0,w,x) A.*exp(-(x-x0).^2./w.^2) + 4;
% myfunc = @(a,b,c,x) a.*(x-b).^2 + c;
% 
% myfit = fittype(@(a,b,c,x) myfunc(a,b,c,x),'independent',{'x'},...
%         'coefficients',{'a','b','c'});
%     opt = fitoptions(myfit);
% opt.StartPoint = [-1.5 60 40];
% 
% % 
%     myfunc = @(x0,y0,a2,a4,x) a2.*(x-x0).^2 + a4.*(x-x0).^4 + y0;    
% myfit = fittype(@(x0,y0,a2,a4,x) myfunc(x0,y0,a2,a4,x),'independent',{'x'},...
%         'coefficients',{'x0','y0','a2','a4'});
% opt = fitoptions(myfit);
% [d0,ind]=min(d);
% f0 = f(ind);
% opt.StartPoint = [d0 f0 1e-3 0];
% 
% fout = fit(f(9:(end-5))',d(9:(end-5))',myfit,opt);
% 
% xx = 56:1:150;
% yy = myfunc(fout.a,fout.b,fout.c,xx);
% 
% x2=[20:1:56];
% y2= -0.0084.*x2+0.9340;

