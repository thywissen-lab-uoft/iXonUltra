figs=get(groot,'Children');
disp(' ');
disp('Closing all non GUI figures.');
for kk=1:length(figs)
   if ~isequal(figs(kk).Tag,'GUI')
       disp(['Closing figure ' num2str(figs(kk).Number) ' ' figs(kk).Name]);
      close(figs(kk)) 
   end
end
disp(' ');


%% Z lattice X Dir 2023/09

runs=[
    2023 09 07 02;
    2023 09 07 03;
    2023 09 07 04;
    2023 09 07 05;
    2023 09 06 08;
    ];
depths = [1.3 3.1 3.1 4 2.2];
data_label = 'Z Lattice XDir + XDT'; 
dVar = 'Yc';
varname = 'zlattice_x';
n = 100;
Tguess = [30 20 20 20 20];


%% Z lattice 2023/09 Y

runs=[
    2023 09 07 06;
    2023 09 07 07;
    2023 09 08 02;
    2023 09 08 03;
    ];

Tguess = [30 20 20 20];

depths = [4 3.1 2.2 1.3 ];
data_label = 'Z Lattice YDir + XDT'; 
dVar = 'Xc';
varname = 'zlattice_y';
n = 110;

%% X lattice 2023/09
runs=[
    2023 09 08 04;
    2023 09 08 05;
    2023 09 08 06;
    2023 09 11 7;

    ];

Tguess = [30 23 20 20];
depths = [1.3 2.2 3.1 3.8];
data_label = 'X Lattice + XDT';
dVar = 'Yc';
varname = 'xlattice';
n = 120;

%% Y lattice 2023/09
runs=[
    2023 09 11 02;
    2023 09 11 03;
    2023 09 11 04;
    2023 09 11 05;
    2023 09 11 06;
    ];
Tguess = [20 18 15 14 13];
depths = [1.3 2.2 3.2 4.2 5.1];
data_label = 'Y Lattice + XDT';
dVar = 'Xc';
varname = 'ylattice';
n = 130;

%% Select the data
% runs = runs_x;
% depths = depths_zx;

%% Load the data
clear data
fname = 'ixon_gaussdata';
[all_data,dirNames,dirDates] = ixon_loadBulk(runs,[fname '.mat']);
data = [all_data.(fname)];
%% 

clear hFs
cmaps = hsv(length(data));
nPlotMax = 6;

clear hFs
j=1;

clear f
clear f_err

out = struct;

for nn=1:length(data)      
    myco = cmaps(nn,:);

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
%         co=get(gca,'colororder');
        t=uicontrol('style','text','string',data_label,'units','pixels',...
            'backgroundcolor','w','horizontalalignment','left','fontsize',10);
        t.Position(3:4)=[hFs(j).Position(3) t.Extent(4)];
        t.Position(1:2)=[5 hFs(j).Position(4)-t.Position(4)];
%         resizeFig(hFs(j),t);
        j=j+1;
    end    
    
     
    % Make Axis and Plot Data
    subplot(3,2,mod(nn-1,nPlotMax)+1);
    plot(X,Y,'o','markerfacecolor',myco,...
        'markeredgecolor',myco*.5,'color',myco,...
        'linewidth',1,'markersize',6);    
    ylabel(dVar);
    xlabel(data(nn).xVar,'interpreter','none');
    

    % Guess the amplitude and offset
    gA=0.5*range(Y);
    gD=(max(Y)+min(Y))*.5;
    
    % Guess the period
%     gB=18;
%     gC=0;
    
%     gB=21-1*nn;
     gC=-2;
    gB = Tguess(nn);
    gE = range(X);

    cosFit=fittype('A*cos(2*pi*t/B+C)*exp(t/E)+D','independent',{'t'},...
        'coefficients',{'A','B','C','D','E'});
    options=fitoptions(cosFit);          
    set(options, 'TolFun', 1E-14);
    set(options, 'StartPoint', [gA, gB, gC,gD, gE]);     
    set(options, 'MaxIter',5000);
    set(options, 'MaxFunEvals',5000);
    set(options,'TolFun',10^-9);
    
    fout = fit(X,Y,cosFit,options);
    cint = confint(fout,0.66);
    hold on
    tt = linspace(0,max(X),100);
    plot(tt,feval(fout,tt),'k-');
    
    str = [num2str(depths(nn)) ' Er, ' num2str(1e3/fout.B,4) ' Hz'];
    f(nn) = 1e3/fout.B;
    f_err(nn) = abs(1e3/cint(2,2)-f(nn));    
    cint(2,2)
     text(.98,.98,str,'units','normalized','verticalalignment','cap',...
        'horizontalalignment','right','interpreter','latex','fontsize',12);   
    out.depths = depths;
    out.freq = f;
    out.freq_err = f_err;
    out.Label = data_label;
    
    assignin('base',varname,out)

end


figure(200);
errorbar(depths,f,f_err,'o','markerfacecolor',[.5 .5 .5],...
    'markeredgecolor','k','color','k',...
    'linewidth',2,'markersize',8); 
xlabel('lattice depth (Er)');
ylabel('measured oscillation frequency (Hz)');
set(gca,'fontsize',14,'xgrid','on','ygrid','on');
xlim([0 6]);
ylim([0 100]);

%%

figure(201);
clf
co = get(gca,'colororder');

uy = [ylattice.depths];
fy = [ylattice.freq];
fey = [ylattice.freq_err];

ux = [xlattice.depths];
fx = [xlattice.freq];
fex = [xlattice.freq_err];

uzx = [zlattice_x.depths];
fzx = [zlattice_x.freq];
fzex = [zlattice_x.freq_err];


uzy = [zlattice_y.depths];
fzy = [zlattice_y.freq];
fzey = [zlattice_y.freq_err];

errorbar([ylattice.depths],[ylattice.freq],[ylattice.freq_err],'o','markerfacecolor',co(1,:),...
    'markeredgecolor',co(1,:)*.5,'color',co(1,:)*.5,...
    'linewidth',1,'markersize',8); 
hold on

errorbar([xlattice.depths],[xlattice.freq],[xlattice.freq_err],'o','markerfacecolor',co(2,:),...
    'markeredgecolor',co(2,:)*.5,'color',co(2,:)*.5,...
    'linewidth',1,'markersize',8); 

errorbar([zlattice_x.depths],[zlattice_x.freq],[zlattice_x.freq_err],'o','markerfacecolor',co(3,:),...
    'markeredgecolor',co(3,:)*.5,'color',co(3,:)*.5,...
    'linewidth',1,'markersize',8); 

errorbar([zlattice_y.depths],[zlattice_y.freq],[zlattice_y.freq_err],'o','markerfacecolor',co(4,:),...
    'markeredgecolor',co(4,:)*.5,'color',co(4,:)*.5,...
    'linewidth',1,'markersize',8); 

legend({'y+xdt','x+xdt','z+xdt x dir','z+xdt y dir'});

xlabel('lattice depth (Er)');
set(gca,'fontsize',14,'xgrid','on','ygrid','on');
xlim([0 6]);
ylim([0 100]);
ylabel('measured oscillation frequency (Hz)');


figure(202);
clf
co = get(gca,'colororder');

fxdt = 33;
errorbar([ylattice.depths],sqrt([ylattice.freq].^2-fxdt^2),[ylattice.freq_err],'o','markerfacecolor',co(1,:),...
    'markeredgecolor',co(1,:)*.5,'color',co(1,:)*.5,...
    'linewidth',1,'markersize',8); 
hold on

errorbar([xlattice.depths],sqrt([xlattice.freq].^2-fxdt^2),[xlattice.freq_err],'o','markerfacecolor',co(2,:),...
    'markeredgecolor',co(2,:)*.5,'color',co(2,:)*.5,...
    'linewidth',1,'markersize',8); 

errorbar([zlattice_x.depths],sqrt([zlattice_x.freq].^2-fxdt^2),[zlattice_x.freq_err],'o','markerfacecolor',co(3,:),...
    'markeredgecolor',co(3,:)*.5,'color',co(3,:)*.5,...
    'linewidth',1,'markersize',8); 

errorbar([zlattice_y.depths],sqrt([zlattice_y.freq].^2-fxdt^2),[zlattice_y.freq_err],'o','markerfacecolor',co(4,:),...
    'markeredgecolor',co(4,:)*.5,'color',co(4,:)*.5,...
    'linewidth',1,'markersize',8); 

legend({'y+xdt','x+xdt','z+xdt x dir','z+xdt y dir'});

xlabel('lattice depth (Er)');
ylabel('inferred trap frequency (Hz)');

set(gca,'fontsize',14,'xgrid','on','ygrid','on');
xlim([0 6]);
ylim([0 100]);



save('trap_freq','uy','fy','fey','ux','fx','fex','uzx','fzx','fzex','uzy','fzy','fzey');
%% Plot and Analyze

% 
% %% UPload data
% doUpload = 1;
% 
% GDrive_root = 'G:\My Drive\Lattice Shared\SharedData\Composite P-wave';
% 
% if  doUpload && exist(GDrive_root,'dir')   
%     gFile = [GDrive_root filesep out_name]; 
%     save(gFile,'data_process','data_out');
%     saveas(hf1,[GDrive_root filesep data_label '_shifts.png'])
%     saveas(hf2,[GDrive_root filesep data_label '_shapes.png'])
%     
%     for jj=1:length(hFs)
%         saveas(hFs(jj),[GDrive_root filesep hFs(jj).Name '.png'])
%     end
% 
% end
