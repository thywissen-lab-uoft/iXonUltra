%% Load in the data

for i=1:length(qgmdata_stripes)

    L1(i) = qgmdata_stripes(i).LatticeDig(1);
    L2(i) = qgmdata_stripes(i).LatticeDig(2);
    L3(i) = qgmdata_stripes(i).LatticeDig(3);
    L4(i) = qgmdata_stripes(i).LatticeDig(4);
    L5(i) = qgmdata_stripes(i).LatticeDig(5);
    L6(i) = qgmdata_stripes(i).LatticeDig(6);

    xC1(i) = L1(i).Xc;
    xC2(i) = L2(i).Xc;
    xC3(i) = L3(i).Xc;
    xC4(i) = L4(i).Xc;
    xC5(i) = L5(i).Xc;
    xC6(i) = L6(i).Xc;

    N1(i) = L1(i).Natoms;
    N2(i) = L2(i).Natoms;
    N3(i) = L3(i).Natoms;
    N4(i) = L4(i).Natoms;
    N5(i) = L5(i).Natoms;
    N6(i) = L6(i).Natoms;

    N2_5(i) = sum([N2(i),N3(i),N4(i),N5(i)]);


end

average =  mean([xC2;xC3;xC4;xC5],1);
average_w = sum([N2.*xC2;N3.*xC3;N4.*xC4;N5.*xC5])./N2_5;

%% Fit all of the stripes

fitf = 'a*sin(2*pi*50*(x+150)/1000-b) + c';
f2 = fit(xvals',xC2',fitf,'Start', [5 pi 115])
f3 = fit(xvals',xC3',fitf,'Start', [5 pi 115])
f4 = fit(xvals',xC4',fitf,'Start', [5 pi 115])
f5 = fit(xvals',xC5',fitf,'Start', [5 pi 115])

ffull = fit(xvals',average',fitf,'Start', [5 pi 115])
ffull_w = fit(xvals',average_w',fitf,'Start',[5 pi 115])

ffull_img = fit(xvals',digdata.Xc',fitf,'Start',[5 pi 115])

%% Plot with everything
hFme = figure(2000);

hFme.Color='w';


xplot = [50:0.1:110];

plot(xvals,average,'oc')

hold on;

plot(xvals,average_w,'ok')

plot(xvals,digdata.Xc,'om')

% plot(xvals,xC1,'o')
plot(xvals,xC2,'oy')
plot(xvals,xC3,'og')
plot(xvals,xC4,'or')
plot(xvals,xC5,'ob')

plot(xplot,f2(xplot),'y')
plot(xplot,f3(xplot),'g')
plot(xplot,f4(xplot),'r')
plot(xplot,f5(xplot),'b')
plot(xplot,ffull(xplot),'c')
plot(xplot,ffull_w(xplot),'k')
plot(xplot,ffull_img(xplot),'m')

xlabel(string(digdata.xVar),Interpreter="none")
ylabel('xC (sites)')

legend('Stripe Mean','Weighted Mean','Full Image','Stripe 2','Stripe 3','Stripe 4','Stripe 5')

hold off;

% filename=fullfile(ixon_imgdir,'figures','digdata1.mat');

%% Make arrays of fit parameters

a = zeros(1,3);
b = zeros(1,3);
c = zeros(1,3);

a_err = zeros(1,3);
b_err = zeros(1,3);
c_err = zeros(1,3);

for i=2:5

    a(i-1) = eval('f'+string(i)+'.a');
    b(i-1) = eval('f'+string(i)+'.b');
    c(i-1) = eval('f'+string(i)+'.c');

    ERR = confint(eval('f'+string(i)));

    a_err(i-1) = mean(abs(ERR(:,1)-a(i-1)));
    b_err(i-1) = mean(abs(ERR(:,2)-b(i-1)));
    c_err(i-1) = mean(abs(ERR(:,3)-c(i-1)));
end

ERR = confint(ffull);
a(end+1) = ffull.a;
b(end+1) = ffull.b;
c(end+1) = ffull.c;

a_err(end+1) = mean(abs(ERR(:,1)-a(end)));
b_err(end+1) = mean(abs(ERR(:,2)-b(end)));
c_err(end+1) = mean(abs(ERR(:,3)-c(end)));

ERR = confint(ffull_w);
a(end+1) = ffull_w.a;
b(end+1) = ffull_w.b;
c(end+1) = ffull_w.c;

a_err(end+1) = mean(abs(ERR(:,1)-a(end)));
b_err(end+1) = mean(abs(ERR(:,2)-b(end)));
c_err(end+1) = mean(abs(ERR(:,3)-c(end)));

ERR = confint(ffull_img);
a(end+1) = ffull_img.a;
b(end+1) = ffull_img.b;
c(end+1) = ffull_img.c;

a_err(end+1) = mean(abs(ERR(:,1)-a(end)));
b_err(end+1) = mean(abs(ERR(:,2)-b(end)));
c_err(end+1) = mean(abs(ERR(:,3)-c(end)));

%% Plot fit parameters

fitfig = figure(2002);
fitfig.Position = [100 100 600 600];
fitfig.Color = 'w';
subplot(4,1,1)

errorbar([2:8],a,a_err)
xlim([1.5 8.5])
xticklabels({'Stripe 2','Stripe 3','Stripe 4','Stripe 5','Stripe Average','Weighted Stripe Average','Full Image'})

ylabel('Oscillation amplitude (sites)')

subplot(4,1,2)

errorbar([2:8],b/(2*pi),b_err/(2*pi))
xlim([1.5 8.5])
xticklabels({'Stripe 2','Stripe 3','Stripe 4','Stripe 5','Stripe Average','Weighted Stripe Average','Full Image'})

ylabel('Phase (2pi)')

subplot(4,1,3)

errorbar([2:8],c,c_err)
xlim([1.5 8.5])
xticklabels({'Stripe 2','Stripe 3','Stripe 4','Stripe 5','Stripe Average','Weighted Stripe Average','Full Image'})

ylabel('COM offset (sites)')

subplot(4,1,4)

errorbar([2:5],mean([N2;N3;N4;N5]./N2_5,2),std([N2;N3;N4;N5]./N2_5,0,2))
xlim([1.75 5.25])
ylim([0 0.4])
xticks([2:5])
xticklabels({'Stripe 2','Stripe 3','Stripe 4','Stripe 5'})

ylabel('Stripe Weight')




%% Plot each fit separately

fplot = figure(2005);
fplot.Color='w';

subplot(4,2,1)
plot(xvals,xC2,'oy')
hold on;
plot(xplot,f2(xplot),'y')
xlim([min(xvals) max(xvals)])
xlabel(string(digdata.xVar),Interpreter="none")
ylabel('xC (sites)')
title('Stripe 2')


subplot(4,2,2)
plot(xvals,xC3,'og')
hold on;
plot(xplot,f3(xplot),'g')
xlim([min(xvals) max(xvals)])
xlabel(string(digdata.xVar),Interpreter="none")
ylabel('xC (sites)')
title('Stripe 3')

subplot(4,2,3)
plot(xvals,xC4,'or')
hold on;
plot(xplot,f4(xplot),'r')
xlim([min(xvals) max(xvals)])
xlabel(string(digdata.xVar),Interpreter="none")
ylabel('xC (sites)')
title('Stripe 4')

subplot(4,2,4)
plot(xvals,xC5,'ob')
hold on;
plot(xplot,f5(xplot),'b')
xlim([min(xvals) max(xvals)])
xlabel(string(digdata.xVar),Interpreter="none")
ylabel('xC (sites)')
title('Stripe 5')

subplot(4,2,5)
plot(xvals,average,'oc')
hold on;
plot(xplot,ffull(xplot),'c')
xlim([min(xvals) max(xvals)])
xlabel(string(digdata.xVar),Interpreter="none")
ylabel('xC (sites)')
title('Stripe Average')

subplot(4,2,6)
plot(xvals,average_w,'ok')
hold on;
plot(xplot,ffull_w(xplot),'k')
xlim([min(xvals) max(xvals)])
xlabel(string(digdata.xVar),Interpreter="none")
ylabel('xC (sites)')
title('Weighted Stripe Average')


subplot(4,2,7)
plot(xvals,digdata.Xc,'om')
hold on;
plot(xplot,ffull_img(xplot),'m')
xlim([min(xvals) max(xvals)])
xlabel(string(digdata.xVar),Interpreter="none")
ylabel('xC (sites)')
title('Full Image')


