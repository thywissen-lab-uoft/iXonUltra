% figs=get(groot,'Children');
% disp(' ');
% disp('Closing all non GUI figures.');
% for kk=1:length(figs)
%    if ~isequal(figs(kk).Tag,'GUI')
%        disp(['Closing figure ' num2str(figs(kk).Number) ' ' figs(kk).Name]);
%       close(figs(kk)) 
%    end
% end
% disp(' ');
%% Quench Data 1Er
% 
% runs=[
%     2023 08 28 03;
%     2023 08 28 04;
%     2023 08 28 05;
%     2023 08 28 06;
%     2023 08 28 07;
%     2023 08 28 08;
%     2023 08 28 09;
%     2023 08 28 10;
%     2023 08 28 11;
%     2023 08 28 12;
%     2023 08 28 13;
%     2023 08 28 14;
%     2023 08 28 15;
%     2023 08 28 16;
%     ];
% 
% fit_types = {'cos','cos','cos','cos','cos','cos','cos','cos','cos','cos','cos','cos','cos','exp','exp','exp','exp','exp'};
% dVar = 'Xc';
% varname = 'quench_1Er_4v4v_08_28';
% n = 100;
% 
% 
% clear data
% fname = 'ixon_gaussdata';
% [all_data,dirNames,dirDates] = ixon_loadBulk(runs,[fname '.mat']);
% data = [all_data.(fname)];
% 
% for kk=1:length(data)
%    data(kk).X =  data(kk).X -5;
%    inds = data(kk).X>=0;
%    data(kk).X = data(kk).X(inds);
%   data(kk).Xc = data(kk).Xc(inds);
% 
% end
%% Quench Data 1
% 
% runs=[
%     2023 09 12 16;
%     2023 09 12 17;
%     2023 09 12 06;
%     2023 09 12 07;
%     2023 09 12 08;
%     2023 09 12 09;
%     2023 09 12 10;
%     2023 09 12 11;
%     2023 09 12 12;
%     2023 09 12 13;
%     2023 09 12 14;
%     2023 09 12 15;
%     2023 09 12 18;
%     2023 09 12 19;
%     2023 09 12 20;
%     2023 09 12 21;
%     2023 09 12 22;
%     2023 09 12 23;
% 
%     ];
% 
% fit_types = {'cos','cos','cos','cos','cos','cos','cos','cos','cos','cos','cos','exp','cos','cos','cos','cos','exp','exp'};
% data_label = '1 Er Quench'; 
% dVar = 'Xc';
% varname = 'quench_1Er_4v4v_09_12';
% n = 100;
% 
% 
% clear data
% fname = 'ixon_gaussdata';
% [all_data,dirNames,dirDates] = ixon_loadBulk(runs,[fname '.mat']);
% data = [all_data.(fname)];
%% Quench Data 1 Er 2v 2v 09/13
% % 
% runs=[
%     2023 09 13 02;
%     2023 09 13 03;
%     2023 09 13 04;
%     2023 09 13 05;
%     2023 09 13 06;
%     2023 09 13 07;
%     2023 09 13 08;
%     2023 09 13 09;
%     2023 09 13 10;
%     2023 09 13 11;
%     2023 09 13 12;
%     2023 09 13 13;
%     2023 09 13 14;
%     2023 09 13 15;
%     2023 09 13 16;
%     2023 09 13 16;
%     2023 09 13 17;
%     2023 09 13 18;
%     2023 09 13 19;
% 
%     ];
% 
% fit_types = {'cos','cos','cos','cos','cos','cos','cos','cos','cos','cos','cos','cos','cos','cos','cos','cos','cos','cos','cos'};
% data_label = '1 Er Quench 2v'; 
% dVar = 'Xc';
% varname = 'quench_1Er_2v2v_09_13';
% n = 100;
% 
% clear data
% fname = 'ixon_gaussdata';
% [all_data,dirNames,dirDates] = ixon_loadBulk(runs,[fname '.mat']);
% data = [all_data.(fname)];
%% Quench Data 3

runs=[
    2023 09 28 08;
    2023 09 28 09;
    2023 09 28 10;
    2023 09 28 11;
    2023 09 28 12;% 
    2023 09 28 13;% 
    2023 09 28 14;% 
    2023 09 28 15;% 
    2023 09 28 16;% 
    2023 09 28 17;% 

    ];

fit_types = {'exp','exp','exp','exp','exp','exp','exp','exp','exp','exp','exp'};
data_label = '2.5 Er Quench'; 
dVar = 'Xc';
varname = 'quench_25Er_4v4v_09_28';
n = 400;
% 
% 
clear data
fname = 'ixon_gaussdata';
[all_data,dirNames,dirDates] = ixon_loadBulk(runs,[fname '.mat']);
data = [all_data.(fname)];
%% Quench Data 3
% 
% runs=[
%     2023 09 28 18;
%     2023 09 28 19;
%     2023 09 28 20;
%     2023 09 28 21;
%     2023 09 28 22;% 
%     2023 09 28 23;% 
%     2023 09 28 24;% 
%     2023 09 28 25;% 
%     2023 09 28 26;% 
%     2023 09 28 27;% 
% 
%     ];
% 
% fit_types = {'cos','cos','cos','cos','cos','cos','cos','cos','exp','exp','exp'};
% data_label = '1.5 Er Quench'; 
% dVar = 'Xc';
% varname = 'quench_15Er_4v4v_09_28';
% n = 500;
% % 
% % 
% clear data
% fname = 'ixon_gaussdata';
% [all_data,dirNames,dirDates] = ixon_loadBulk(runs,[fname '.mat']);
% data = [all_data.(fname)];
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
out = struct;

for nn=1:length(data)      
    myco = cmaps(nn,:);
    N(nn) = mean(data(nn).N);
    N_err(nn) = std(data(nn).N);
    b = data(nn).Params.conductivity_FB_field;
%     b = b-0.11;    
    B(nn) = b;
    a(nn) = B2a(b);    
    % X Data
    X = data(nn).X;
    % Ydata
    Y = data(nn).(dVar); 
    
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
         
    % Make Axis and Plot Data
    subplot(3,3,mod(nn-1,nPlotMax)+1);
    plot(X,Y,'o','markerfacecolor',myco,...
        'markeredgecolor',myco*.5,'color',myco,...
        'linewidth',1,'markersize',4);    
    hold on;
    ylabel(dVar);
    xlabel(data(nn).xVar,'interpreter','none');
    % amplitude
    gA=0.5*range(Y);          
    % offset
    gD=(max(Y)+min(Y))*.5;    
    % phase
    gC = 47;    
    % Period
    gB = 23;    
    % Decay rate
    gE = 1/10;
    
    switch fit_types{nn}
        case 'cos'
%             inds = X>=7.5;
%             X=X(inds);
%             Y=Y(inds);
            
            cosFit=fittype('A*cos(2*pi*t/B+C)*exp(-E/2*t)+D','independent',{'t'},...
                'coefficients',{'A','B','C','D','E'});
            options=fitoptions(cosFit);          
            set(options, 'TolFun', 1E-14);
            set(options, 'StartPoint', [gA, gB, gC,gD, gE]);     
            set(options, 'MaxIter',10000);
            set(options, 'MaxFunEvals',10000);
            set(options,'Robust','Bisquare');
            set(options,'TolFun',10^-9);
            options.Upper = [40 25 500 270 200];  
            fout = fit(X,Y,cosFit,options);
            cint = confint(fout,0.66);
            hold on
            tt = linspace(0,max(X),100);
            plot(tt,feval(fout,tt),'k-');
            f(nn) = 1e3/fout.B;
            f_err(nn) = abs(1e3/cint(2,2)-f(nn));        
            gamma(nn) = fout.E;     
            gamma_err(nn) = range(cint(:,5))/2;
            tau(nn) = 2/fout.E;
            
        case 'exp'
            tt = linspace(0,max(X),150);

            expFit=fittype(@(A,E,D,t) A*exp(-E/2*t)+D,'independent',{'t'},...
                'coefficients',{'A','E','D'});
            options=fitoptions(expFit);          
            set(options, 'TolFun', 1E-14);
            set(options, 'StartPoint', [-gA 1/10 262]);     
            set(options, 'MaxIter',5000);
            set(options, 'MaxFunEvals',5000);
            set(options,'TolFun',10^-9);
            fout = fit(X,Y,expFit,options);
            plot(tt,feval(fout,tt),'k-');
            f(nn) = NaN;
            f_err(nn) = NaN;
            gamma(nn) = fout.E;
            gamma_err(nn) = range(cint(:,5))/2;
            tau(nn) = 2/fout.E;
        case 'oexp'
            tt = linspace(0,max(X),150);

            func = @(A,B,G,omega,D,t) ...
                A*exp(t*(-G/2+sqrt((G/2)^2-omega^2))) + ...
                B*exp(t*(-G/2-sqrt((G/2)^2-omega^2))) + ...
                D;

            expFit=fittype(@(A,G,D,t) func(A,0,G,2*pi*0.04,D,t),'independent',{'t'},...
                'coefficients',{'A','G','D'});
            options=fitoptions(expFit);          
            set(options, 'TolFun', 1E-14);
            set(options, 'StartPoint', [-gA 1 270]);     
            set(options, 'MaxIter',5000);
            set(options, 'MaxFunEvals',5000);
            set(options,'TolFun',10^-9);
            fout = fit(X,Y,expFit,options);
            plot(tt,feval(fout,tt),'k-');
            f(nn) = NaN;
            f_err(nn) = NaN;
            gamma(nn) = fout.G;
            gamma_err(nn) = range(cint(:,5))/2;
            tau(nn) = 2/fout.G;
    end
    

     
    str = [num2str(b) ' G,' '$\tau =' num2str(tau(nn),3) '$ms'];
    text(.98,.98,str,'units','normalized','verticalalignment','cap',...
         'horizontalalignment','right','interpreter','latex','fontsize',12,...
         'interpreter','latex');       
    str2 = '$\exp(-\Gamma/2 t),\tau = 2/\Gamma$';
    text(.02,.02,str2,'units','normalized','verticalalignment','bottom',...
         'horizontalalignment','left','interpreter','latex','fontsize',12,...
         'interpreter','latex');     
end
%% Summarize the data

out.field = B;
out.freq = f;
out.freq_err = f_err;
out.gamma = gamma;
out.gamma_err = gamma_err;
out.tau = tau;
out.N = N;
out.N_err = N_err;
% out.Label = varname;    

field = B;
freq = f;
freq_err = f_err;
gamma = gamma;
gamma_err= gamma_err;
tau = tau;
N = N;  
assignin('base',varname,out) 

%% 

% Counts Track
hf1=figure(1001);
set(gcf,'WindowStyle','docked','color','w');
clf
p = get(gcf,'Position');
errorbar(B,N,N_err,'o','markerfacecolor',[.5 .5 .5],...
    'markeredgecolor','k','color','k',...
    'linewidth',2,'markersize',8); 
xlabel('field (G)');
ylabel('fluoresence counts (arb)');
set(gca,'fontsize',14,'xgrid','on','ygrid','on');
t=uicontrol('style','text','string',varname,'units','pixels',...
        'backgroundcolor','w','horizontalalignment','left','fontsize',10);
t.Position(3:4)=[p(3) t.Extent(4)];
t.Position(1:2)=[5 5];
yL = get(gca,'Ylim');
ylim([0 yL(2)]);
% ylim([0 6e7]);

% Frequency Track
hf2=figure(1002);
set(gcf,'WindowStyle','docked','color','w');
clf
p = get(gcf,'Position');
errorbar(B,f,f_err,'o','markerfacecolor',[.5 .5 .5],...
    'markeredgecolor','k','color','k',...
    'linewidth',2,'markersize',8); 
xlabel('field (G)');
ylabel('measured oscillation frequency (Hz)');
set(gca,'fontsize',14,'xgrid','on','ygrid','on'); 
t=uicontrol('style','text','string',varname,'units','pixels',...
        'backgroundcolor','w','horizontalalignment','left','fontsize',10);
t.Position(3:4)=[p(3) t.Extent(4)];
t.Position(1:2)=[5 p(4)-t.Position(4)];
    
% rate vs field Track
hf3=figure(1003);
set(gcf,'WindowStyle','docked','color','w');
clf
p = get(gcf,'Position');
t=uicontrol('style','text','string',varname,'units','pixels',...
'backgroundcolor','w','horizontalalignment','left','fontsize',10);
t.Position(3:4)=[p(3) t.Extent(4)];
t.Position(1:2)=[5 p(4)-t.Position(4)];
errorbar(B,gamma,gamma_err,'o','markerfacecolor',[.5 .5 .5],...
     'markeredgecolor','k','color','k',...
     'linewidth',2,'markersize',8); 
xlabel('field (G)');
ylabel('decay rate \Gamma (ms^{-1})');
set(gca,'fontsize',12,'xgrid','on','ygrid','on');
 yl=get(gca,'ylim');
 ylim([0 yl(2)]);
 
% rate vs field Track
hf4=figure(1004);
set(gcf,'WindowStyle','docked','color','w');
clf
p = get(gcf,'Position');
t=uicontrol('style','text','string',varname,'units','pixels',...
'backgroundcolor','w','horizontalalignment','left','fontsize',10);
t.Position(3:4)=[p(3) t.Extent(4)];
t.Position(1:2)=[5 p(4)-t.Position(4)];
 errorbar(a.^2,gamma,gamma_err,'o','markerfacecolor',[.5 .5 .5],...
     'markeredgecolor','k','color','k',...
     'linewidth',2,'markersize',8); 
xlabel('a^2 (a_0)');
ylabel('decay rate \Gamma (ms^{-1})');
set(gca,'fontsize',12,'xgrid','on','ygrid','on');
 yl=get(gca,'ylim');
 ylim([0 yl(2)]);
 xl=get(gca,'xlim');
 xlim([0 xl(2)]);
%% Save and UPload data
doUpload = 0;
GDrive_root =['G:\My Drive\Lattice Shared\SharedData\Conductivity_2023' filesep varname];
doSave = 1;
 
 if doSave 
    curpath = fileparts(mfilename('fullpath'));
    save(fullfile(curpath,varname),'field','freq','freq_err','gamma','gamma_err','tau','N','N_err');
 end
 
 try
     mkdir(GDrive_root);
 end
 
 if  doUpload && exist(GDrive_root,'dir')   
     disp('upload to google drive');
    gFile = [GDrive_root filesep varname '.mat']; 
    save(gFile,varname);
    saveas(hf1,[GDrive_root filesep varname '_counts.png'])
    saveas(hf2,[GDrive_root filesep varname '_freq.png'])
    saveas(hf3,[GDrive_root filesep varname '_decayB.png'])
    saveas(hf4,[GDrive_root filesep varname '_decaya2.png'])

    for jj=1:length(hFs)
        saveas(hFs(jj),[GDrive_root filesep varname '_' hFs(jj).Name '.png'])
    end
 
 end
