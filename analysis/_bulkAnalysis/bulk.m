
runs=[
    2023 10 05 17;
    2023 10 05 18;
    2023 10 05 19;
    2023 10 05 20;
    2023 10 05 21;
    2023 10 05 22;
    2023 10 05 23;
    2023 10 05 24;
    2023 10 05 25;
    2023 10 05 26;

    ];

fit_types = {'exp','exp','exp','exp','exp','exp','exp','exp','exp','exp','exp'};
data_label = '2.5 Er Quench'; 
dVar = 'Xc';
varname = 'shake_25Er_4v4v_70hz';
n = 400;
% 
%%


runs=[
    2023 10 06 02;
    2023 10 06 03;
    2023 10 06 04;
    2023 10 06 05;
    2023 10 06 06;
    2023 10 06 07;
    2023 10 06 08;
    2023 10 06 09;

    ];

fit_types = {'exp','exp','exp','exp','exp','exp','exp','exp','exp','exp','exp'};
data_label = '2.5 Er Quench'; 
dVar = 'Xc';
varname = 'shake_25Er_4v4v_195G';
n = 400;
% 
% 
clear data
fname = 'ixon_gaussdata';
[all_data,dirNames,dirDates] = ixon_loadBulk(runs,[fname '.mat']);
data = [all_data.(fname)];


%%
runs=[
    2023 10 06 12;
    2023 10 06 13;
    2023 10 06 14;
    2023 10 06 15;

    ];

fit_types = {'exp','exp','exp','exp','exp','exp','exp','exp','exp','exp','exp'};
data_label = '2.5 Er Quench'; 
dVar = 'Xc';
varname = 'shake_25Er_4v4v_20G';
n = 400;
% 
% 
clear data
fname = 'ixon_gaussdata';
[all_data,dirNames,dirDates] = ixon_loadBulk(runs,[fname '.mat']);
data = [all_data.(fname)];
%%
clear data
fname = 'ixon_gaussdata';
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
out = struct;

for nn=1:length(data)      
    myco = cmaps(nn,:);
    N(nn) = mean(data(nn).N);
    N_err(nn) = std(data(nn).N);
    
%     B(nn)= data(nn).Params.conductivity_FB_field_maybe_calibrated;
    B(nn)= 20;

f(nn)= data(nn).Params.conductivity_mod_freq;
    a(nn) = B2a(B(nn));    
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

    title([num2str(B(nn)) ' G, ' num2str(f(nn)) ' Hz']);
        ylim([270 290]);
end
