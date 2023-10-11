runs=[
    2023 10 07 02;
    2023 10 07 03;
    2023 10 07 04;
    2023 10 07 05;
    2023 10 07 06;
    2023 10 07 07;
    2023 10 07 08;
    2023 10 07 09;
    2023 10 07 10;
    2023 10 07 11;
    2023 10 07 12;
    2023 10 07 13;
    2023 10 07 14;
    ];

% fit_types = {''};
data_label = '2.5 Er'; 
dVar = 'Xc';
varname = 'shake_25Er_2v2v_170G_73mW';
n = 400;
% 
% %
clear data
fname = 'digdata';
[all_data,dirNames,dirDates] = ixon_loadBulk(runs,[fname '.mat']);
data = [all_data.(fname)];
pddir = 'X:\LabJackLogs\ODTQPD\2023\2023.10\10.07';

%%
runs=[
    2023 10 07 22;
    2023 10 07 23;
    2023 10 07 24;
    2023 10 07 25;
    2023 10 07 26;
    2023 10 07 27;
    2023 10 07 28;
    2023 10 07 29;
    2023 10 07 30;
    2023 10 07 31;
    2023 10 07 32;
    ];

% fit_types = {''};
data_label = '2.5 Er Quench'; 
dVar = 'Xc';
varname = 'shake_25Er_2v2v_54Hz_70mW';
n = 400;
% 
% 
clear data
fname = 'digdata';
[all_data,dirNames,dirDates] = ixon_loadBulk(runs,[fname '.mat']);
data = [all_data.(fname)];
pddir = 'X:\LabJackLogs\ODTQPD\2023\2023.10\10.07';

%%

runs = [2023 10 07 33;
    2023 10 07 34;
    2023 10 07 35;
    2023 10 07 36;
    2023 10 07 37;
    2023 10 07 38;
    2023 10 07 39;
    2023 10 07 40;
    2023 10 07 41;
    2023 10 07 42;
    2023 10 07 43;
    2023 10 07 44;
    2023 10 07 45;
    2023 10 07 46;
    2023 10 07 47;
    2023 10 07 48;
    2023 10 07 49;
    2023 10 07 50;
    2023 10 07 51;
    2023 10 07 52;
    ];

% fit_types = {''};
data_label = '2.5 Er Quench'; 
dVar = 'Xc';
varname = 'shake_25Er_2v2v_54Hz_80mW';
n = 400;
% 
% 
clear data
fname = 'digdata';
[all_data,dirNames,dirDates] = ixon_loadBulk(runs,[fname '.mat']);
data = [all_data.(fname)];

%%

runs = [2023 10 09 03;
    2023 10 09 04;
2023 10 09 05;
2023 10 09 06;
2023 10 09 07;
2023 10 09 08;
2023 10 09 09;
2023 10 09 10;
2023 10 09 11;
2023 10 09 12;
2023 10 09 13;
2023 10 09 14;
2023 10 09 15;
    ];

% fit_types = {''};
dVar = 'Xc';
varname = 'shake_25Er_2v2v_54Hz_80mW';
data_label = varname; 

n = 400;
% 
% 
clear data
fname = 'digdata';
[all_data,dirNames,dirDates] = ixon_loadBulk(runs,[fname '.mat']);
data = [all_data.(fname)];
%% 
B2a = @(Bfield) 167*(1-(6.910)./(Bfield-202.15));
clear hFs
cmaps = hsv(length(data));
nPlotMax = 9;

clear hFs
j=1;

clear B
clear f
clear f_err
clear a
clear tau
clear N
clear N_err
clear gamma
clear gamma_err
clear myamp
clear myamp_err

clear phi
clear phi_err

out = struct;

for nn=1:length(data)      
    fme = datestr(data(nn).Params(end).ExecutionDateStr,'YYYY-mm-dd_HH-MM-SS');
   
    pdsrc = 'X:\LabJackLogs\ODTQPD';
    pddir = fullfile(pdsrc,fme(1:4),[fme(1:4) '.' fme(6:7)],[fme(6:7) '.' fme(9:10)]);

    
    
    fname = ['ODTQPD_' fme '.mat'];
    qpdfile = fullfile(pddir,fname);
       
    if ~exist(qpdfile)
        warning('cant find pd file');
        qpdfile = fullfile(pddir,'ODTQPD_2023-10-07_19-41-04.mat');
    end
      dpd = load(qpdfile);
   

    myco = cmaps(nn,:);
    N(nn) = mean(data(nn).N);
    N_err(nn) = std(data(nn).N);
    
    B(nn)= data(nn).Params.conductivity_FB_field_maybe_calibrated;
    %B(nn)= 20;

    f(nn)= data(nn).Params.conductivity_mod_freq;
    a(nn) = B2a(B(nn));    
    % X Data
    X = data(nn).X;
    % Ydata
    Y = data(nn).(dVar); 
    N = data(nn).N;
    Xs = data(nn).Xs;
    
    
    binds = isnan(Y);
       X(binds)=[];
     Y(binds)=[];
      N(binds)=[];
      Xs(binds)=[];  
      

    
%      binds = abs(Y-Ym)>8;    
%      X(binds)=[];
%      Y(binds)=[];
%      N(binds)=[];
%      Xs(binds)=[];

      Ym = median(Y);
    binds = N/median(N)<.25;    
     X(binds)=[];
     Y(binds)=[];
     N(binds)=[];
     Xs(binds)=[];
    
%     tpd=1e3*dpd.t-775+X(end);
    tpd = 1e3*dpd.t-680;
    
    vpd1 = dpd.data(:,1)./dpd.data(:,3);
    vpd2 = dpd.data(:,4)./dpd.data(:,6);
    
    % add minus sign to account for the fact that + peizo movies cloud
    % negative
    vpd1 = -vpd1;
    vpd2 = -vpd2;
    
    i1 = find(tpd>=min(X),1);
    i2 = find(tpd>=max(X),1);
    
    t_sub = tpd(i1:i2);
    vpd1_sub = vpd1(i1:i2);
     
    % Make a new figure if necessary
    if ~mod(nn-1,nPlotMax)
        % Plot Data    
        hFs(j)=figure(n+floor(nn/nPlotMax));
        clf
        hFs(j).Color='w';
        hFs(j).Position=[100 50 800 400];
        hFs(j).Name = [dVar '_' num2str(j)];
        hFs(j).WindowStyle='docked';
        t=uicontrol('style','text','string',varname,'units','pixels',...
            'backgroundcolor','w','horizontalalignment','left','fontsize',10);
        t.Position(3:4)=[hFs(j).Position(3) t.Extent(4)];
        t.Position(1:2)=[5 hFs(j).Position(4)-t.Position(4)];
        j=j+1;
    end    
    
    T = 1e3/f(nn);
    

    myfunc = @(A,B,C,t) A*sin(2*pi*t/T + B) + C;    

    myfit = fittype(@(A,B,C,t) myfunc(A,B,C,t),'independent',{'t'},...
        'coefficients',{'A','B','C'});
    
    myfit2 = fittype(@(A,B,C,t) myfunc(A,B,C,t),'independent',{'t'},...
        'coefficients',{'A','B','C'});
    
    Agp = (max(vpd1_sub)-min(vpd1_sub))/2;
    Cgp = mean(vpd1_sub);
    Bgp =  mod(2*pi*(150)/T+pi,2*pi);% additional pi since positive piezo moves cloud to -
    opt = fitoptions(myfit);
    opt.StartPoint = [Agp Bgp Cgp];
    fout_pd = fit(t_sub',vpd1_sub,myfit,opt);
    tt=linspace(min(X),max(X),1e3);
    
    Ag = (max(Y) - min(Y))*.5;
    Cg = mean(Y);
    Bg = fout_pd.B;
    
    opt.StartPoint = [Ag Bg Cg];
    opt.Robust = 'bisquare';
    opt.Lower  = [0 -10*pi Cg-20];
     opt.Upper  = [Ag*5 +10*pi Cg+20];

    fout_dig = fit(X',Y',myfit,opt);
    c=confint(fout_dig,0.667);   
    

         
    % Make Axis and Plot Data
    axbot=subplot(3,3,mod(nn-1,nPlotMax)+1);    
    
    pos = axbot.Position;
    pos1 = pos;
    pos1(4) = pos(4)*.5;
    
    pos2 = pos;
    pos2(2)  = pos1(2)+pos1(4);    
    pos2(4) = pos(4)*.25;  
    
    pos3 = pos;
    pos3(2) = pos2(2)+pos2(4);
    pos3(4) = pos(4)*.25;      
       
    axbot.Position = pos1;        
    axmid = axes('position',pos2);
    axtop = axes('position',pos3);

    % Size Plot
    axes(axtop);
    plot(X,Xs,'^','markerfacecolor',myco,...
        'markeredgecolor',myco*.5,'color',myco,...
        'linewidth',1,'markersize',4);    
    set(gca,'YColor',myco*.5);
    ylabel(['\sigma X (a_L)']);
    xlim([min(X) max(X)]); 
    title([num2str(B(nn)) ' G, ' num2str(f(nn)) ' Hz']);       

    % Atom Number plot
    axes(axmid);
    plot(X,N,'s','markerfacecolor',myco,...
        'markeredgecolor',myco*.5,'color',myco,...
        'linewidth',1,'markersize',4);    
    set(gca,'YColor',myco*.5);
    ylabel(['N']);
    xlim([min(X) max(X)]); 
    

    
    axes(axbot);
%     
    yyaxis right
    plot(tpd,vpd1,'k-');
    hold on
    plot(tt,feval(fout_pd,tt),'-','color',[.5 .5 .5],'linewidth',2);
    set(gca,'YColor','k');

    yyaxis left
    plot(tt,feval(fout_dig,tt),'-','linewidth',2,'color',myco);
    hold on
    plot(X,Y,'o','markerfacecolor',myco,...
        'markeredgecolor',myco*.5,'color',myco,...
        'linewidth',1,'markersize',4);    
    set(gca,'YColor',myco*.5);
    ylabel([dVar ' (a_L)']);
    xlabel(data(nn).xVar,'interpreter','none');
    hold on
    yyaxis left
    plot(tpd,vpd1,'k-');
    hold on
    plot(tt,feval(fout_pd,tt),'-','color',[.5 .5 .5],'linewidth',2);
    set(gca,'YColor','k');   
    ylim(Ym+[-15 15]);
    xlim([min(X) max(X)]); 
    
    
           

    myamp(nn) = fout_dig.A;
    myamp_err(nn) = (c(2,1)-c(1,1))/2;

   % phi(nn) = mod(fout_dig.B,2*pi)-mod(fout_pd.B,2*pi);
    phi(nn) = mod(fout_dig.B-fout_pd.B,2*pi)-1.5*pi;
    %phi(nn) = mod(fout_dig.B-Bgp,2*pi)-pi/2;
    phi_err(nn) = (c(2,2)-c(1,2))/2;

end
%%

clear Nbar
clear Nerr
for nn=1:length(data)
    Nbar(nn) = mean(data(nn).N);
    Nerr(nn) = std(data(nn).N);
end


f1=figure(1213)
clf
f1.WindowStyle='docked';
subplot(131);
errorbar(f,myamp,myamp_err,'ko','markerfacecolor','k');
xlabel('frequency (Hz)');
ylabel('amplitude (lattice sites)');
grid on
ylim([0 5]);
subplot(132);
errorbar(f,phi/(2*pi),phi_err/(2*pi),'o','markerfacecolor','b');
ylabel('phase difference (2*pi)');
ylim([-0.5 0.5]);
grid on
xlabel('frequency (Hz)');

subplot(133);
errorbar(f,Nbar,Nerr,'ko','markerfacecolor','k');
ylabel('N atoms');
xlabel('frequency (Hz)');
grid on

%%

%%
f2=figure(1214)
clf
f2.WindowStyle='docked';

subplot(131);
errorbar(B,myamp,myamp_err,'ko','markerfacecolor','k');
xlabel('Bfield (G)');
ylabel('amplitude (lattice sites)');
grid on
ylim([0 5]);
subplot(132);
errorbar(B,phi/(2*pi),phi_err/(2*pi),'o','markerfacecolor','b');
ylabel('phase difference (2*pi)');
ylim([-0.5 0.5]);
grid on
xlabel('Bfield (G)');

subplot(133);
errorbar(B,Nbar,Nerr,'ko','markerfacecolor','k');
ylabel('N atoms');
xlabel('Bfield (G)');
grid on
%%
f3=figure(125)
clf
f3.WindowStyle='docked';

subplot(131);
errorbar(B2a(B).^2,myamp,myamp_err,'ko','markerfacecolor','k');
xlabel('a^2 (a_0^2)');
ylabel('amplitude (lattice sites)');
grid on
ylim([0 5]);
subplot(132);
errorbar(B2a(B).^2,phi/(2*pi),phi_err/(2*pi),'o','markerfacecolor','b');
ylabel('phase difference (2*pi)');
ylim([-0.5 0.5]);
grid on
xlabel('a^2 (a_0^2)');

subplot(133);
errorbar(B2a(B).^2,Nbar,Nerr,'ko','markerfacecolor','k');
ylabel('N atoms');
xlabel('a^2 (a_0^2)');
grid on
%%

