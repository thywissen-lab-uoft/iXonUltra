

x = [20,30,40,45,50,60,70,80,90,100,110,120];
y = [0.72,0.70,0.71,0.73,0.76,0.87,1.03,1.27,1.56,1.91,2.30,2.73];

% myfunc = @(A,x0,w,x) A.*exp(-(x-x0).^2./w.^2) + 4;
myfunc = @(a,b,c,x) a.*(x-b).^2 + c;

myfit = fittype(@(a,b,c,x) myfunc(a,b,c,x),'independent',{'x'},...
        'coefficients',{'a','b','c'});

opt = fitoptions(myfit);
opt.StartPoint = [-1.5 20 40];
fout = fit(x(1:end)',y(1:end)',myfit,opt);

xx = 20:1:120;
yy = myfunc(fout.a,fout.b,fout.c,xx);

clf(f33)
f33 = figure(33);
plot(x,y,'ko')
hold on;
plot(xx,yy,'r')
% plot([20 40], [1.63,1.63],'r')
% plot([98 100], [4,4],'r')
xlabel('Drive Frequency (Hz)')
ylabel('Drive Amplitude (V)')
