%% 

runs = [
    2025 02 03 02
    2025 02 03 03
    2025 02 03 04
    2025 02 03 05
    2025 02 03 06
    2025 02 04 01
    2025 02 04 02
    2025 02 04 03
    2025 02 04 04
    2025 02 04 05
    2025 02 04 06
    2025 02 04 07
    2025 02 04 08
    2025 02 04 09
    ];
composite_data= struct;
composite_data.Runs = runs;
% composite_data.Name = '2025_01_29 201.1 G 55 Hz 63.7 mW, 1 Er pulse, vary vert disp';
composite_data.Name = '2025_02_04 ODT trap freq. with vertical displacement';
%% Analyze the runs

composite_opts=struct;
qpd_super(composite_data,composite_opts)

%% Analyze the Y-position across different runs

clear Y1all
clear Y2all
clear Y1allErr
clear Y2allErr

hF = figure;
hF.Color='w';
hF.Position=[10 50 1400 1000];
co=parula(length(composite_data));

for super_composite_data_index=1:length(composite_data)
    dir_list = ixon_findRunDirectory(composite_data(super_composite_data_index).Runs);
    
    fprintf('(%.0f/%.0f) Pulling QPD data...',super_composite_data_index,length(composite_data));
    
    clear X;
    clear Y1;
    clear Y2;
    clear Y1Err;
    clear Y2Err;
    
    for super_directory_index=1:length(dir_list)

        filename = fullfile(dir_list{super_directory_index},'figures','digdata.mat');
        bd = load(filename);
%         bd = bd.ixon_boxdata;
        
        % Get the vertical position
%         Vert_disp(super_directory_index) = [bd.Params(1).xdtB_piezo_vert_disp_amplitude];

        X(super_directory_index) = bd.Params(1).ExecutionDate;
        
        % Load the QPD data
        q = load(fullfile(dir_list{super_directory_index},'figures','qpd.mat'));
        qpd_odt = q.QPD_Modulation;
        Y1(super_directory_index) = 1e3*mean([qpd_odt.Y1],"all");
        Y2(super_directory_index) = 1e3*mean([qpd_odt.Y2],"all");
        Y1Err(super_directory_index) = 1e3*std2([qpd_odt.Y1]);
        Y2Err(super_directory_index) = 1e3*std2([qpd_odt.Y2]);
        
    end
    
    fprintf('done! \n')
    
    subplot(2,2,1)
    set(gca,'fontsize',10,'box','on','linewidth',1);
    errorbar(X,Y1,Y1Err,'o','MarkerFaceColor',co(super_composite_data_index,:),...
        'Color',co(super_composite_data_index,:)*.5,'linewidth',1,'markersize',4);
    hold on
    
    subplot(2,2,2)
    set(gca,'fontsize',10,'box','on','linewidth',1);
    errorbar(X,Y2,Y2Err,'o','MarkerFaceColor',co(super_composite_data_index,:),...
        'Color',co(super_composite_data_index,:)*.5,'linewidth',1,'markersize',4);
    hold on
    
    subplot(2,2,3)
    set(gca,'fontsize',10,'box','on','linewidth',1);
    errorbar(Y1,Y2,Y2Err,Y2Err,Y1Err,Y1Err,'o','MarkerFaceColor',co(super_composite_data_index,:),...
        'Color',co(super_composite_data_index,:)*.5,'linewidth',1,'markersize',4);
    hold on
    
    
    Y1all(super_composite_data_index)    = mean(Y1);
    Y2all(super_composite_data_index)    = mean(Y2);
    Y1allErr(super_composite_data_index) = std(Y1);
    Y2allErr(super_composite_data_index) = std(Y2);
    
end

%% Add calibration dates and labels

Xvar = 'ExecutionDate';
subplot(2,2,1)
set(gca,'fontsize',10,'box','on','linewidth',1);
xlabel(Xvar)
ylabel('xdt1 y-position (mV/V)')
xlim([datenum(2024,10,01) datenum(2025,03,01)])
ylim([-180 -120])

% Calibration dates
plot(datenum(2024,10,18)*ones(100),linspace(-180,-120),'r:','linewidth',1)
plot(datenum(2024,12,18)*ones(100),linspace(-180,-120),'r:','linewidth',1)
plot(datenum(2025,01,08)*ones(100),linspace(-180,-120),'r:','linewidth',1)
plot(datenum(2025,01,21)*ones(100),linspace(-180,-120),'r:','linewidth',1)
plot(datenum(2025,01,29)*ones(100),linspace(-180,-120),'r:','linewidth',1)
plot(datenum(2025,02,06)*ones(100),linspace(-180,-120),'r:','linewidth',1)
plot(datenum(2025,02,07)*ones(100),linspace(-180,-120),'r:','linewidth',1)

datetick x

subplot(2,2,2)
set(gca,'fontsize',10,'box','on','linewidth',1);
xlabel(Xvar)
ylabel('xdt2 y-position (mV/V)')
datetick x
xlim([datenum(2024,10,01) datenum(2025,03,01)])
ylim([-70 -60])

%Calibration dates
plot(datenum(2024,10,18)*ones(100),linspace(-70,-60),'r:','linewidth',1)
plot(datenum(2024,12,18)*ones(100),linspace(-70,-60),'r:','linewidth',1)
plot(datenum(2025,01,08)*ones(100),linspace(-70,-60),'r:','linewidth',1)
plot(datenum(2025,01,21)*ones(100),linspace(-70,-60),'r:','linewidth',1)
plot(datenum(2025,01,29)*ones(100),linspace(-70,-60),'r:','linewidth',1)
plot(datenum(2025,02,06)*ones(100),linspace(-70,-60),'r:','linewidth',1)
plot(datenum(2025,02,07)*ones(100),linspace(-70,-60),'r:','linewidth',1)

datetick x

subplot(2,2,3)
set(gca,'fontsize',10,'box','on','linewidth',1);
xlabel('xdt1 y-position (mV/V)')
ylabel('xdt2 y-position (mV/V)')
xlim([-180 -120])
ylim([-70 -60])

%% Plot sum rule ratio

% Load the data
src_Dec05 = "2024_12_05_2.5Er_201.1G_full_spectrum_SHIFT_correct_disp_RESCALE.mat";
src_Dec09 = "2024_12_09_2.5Er_201.1G_full_spectrum_SHIFT_correct_disp_RESCALE.mat";
src_12_16_spectrum = "2024_12_13_2.5Er_full_spectrum_vary_U_SHIFT_correct_disp_RESCALE.mat";
src_Nov11_spectrum = "2024_11_2.5Er_190G_full_spectrum_SHIFT_correct_disp_RESCALE.mat";
src_Nov_vT = "2024_11_2.5Er_201.1G_full_spectrum_vary_temp_SHIFT_correct_disp_RESCALE.mat";
src_Jan_vT = "2025_01_23_2.5Er_201.1G_full_spectrum_vary_temp_SHIFT_correct_disp_RESCALE.mat";
src_Jan27 = "2025_01_27_2.5Er_full_spectrum_vary_U_SHIFT_correct_disp_RESCALE.mat";

src = [src_Nov_vT,src_Nov11_spectrum,src_Dec05,src_Dec09,src_12_16_spectrum,...
        src_Jan_vT,src_Jan27];
    
    
RescaleFactor = [];

for jj=1:length(src)
    s = load(src(jj));
    s = s.spectral_summary_rescaled;
    
    for kk=1:length(s)
        RescaleFactor(end+1) = s(kk).RescaleFactor;
    end
end


% Plot it
subplot(2,2,4)
scatter(Y1all,Y2all,40,RescaleFactor,"filled");
set(gca,'fontsize',10,'box','on','linewidth',1);
xlabel('xdt1 y-position (mV/V)')
ylabel('xdt2 y-position (mV/V)')
c=colorbar();
c.Label.String='Rescale Factor';
xlim([-180 -120])
ylim([-70 -60])



%% Prepare the data

%%%%% Choose what you want to plot %%%%

X = Vert_disp;
Xvar = 'Vertical displacement (V)';

%%%% Other analysis %%%%

% Fit Y to a line
fitY = 1;
if fitY
    Xplot = linspace(min(X),max(X));
    fout1 = fit(X',Y1','poly1');
    c1 = confint(fout1,0.667); 
    legStr1 = sprintf('y = (%.2f ± %.2f) x + (%.2f ± %.2f)',...
        fout1.p1,0.5*(c1(2,1)-c1(1,1)),fout1.p2,0.5*(c1(2,2)-c1(1,2)));

    fout2 = fit(X',Y2','poly1');
    c2 = confint(fout2,0.667); 
    legStr2 = sprintf('y = (%.2f ± %.2f) x + (%.2f ± %.2f)',...
        fout2.p1,0.5*(c2(2,1)-c2(1,1)),fout2.p2,0.5*(c2(2,2)-c2(1,2)));
end

%% Plot the results

hF = figure;
    hF.Color='w';
    hF.Position=[10 50 1400 400];
    co=get(gca,'colororder');
    
    subplot(1,3,1)
    set(gca,'fontsize',10,'box','on','linewidth',1);
    errorbar(X,Y1,Y1Err,'o','MarkerFaceColor',co(1,:),...
        'Color',co(1,:)*.5,'linewidth',1,'markersize',8);
    xlabel(Xvar)
    ylabel('xdt1 y-position (mV/V)')
    if fitY
        hold on
        plot(Xplot,feval(fout1,Xplot),'linewidth',2);
        legend('',legStr1)
    end
    
    subplot(1,3,2)
    set(gca,'fontsize',10,'box','on','linewidth',1);
    errorbar(X,Y2,Y2Err,'o','MarkerFaceColor',co(1,:),...
        'Color',co(1,:)*.5,'linewidth',1,'markersize',8);
    xlabel(Xvar)
    ylabel('xdt2 y-position (mV/V)') 
    if fitY
        hold on
        plot(Xplot,feval(fout2,Xplot),'linewidth',2);
        legend('',legStr2)
    end


    
    subplot(1,3,3)
    set(gca,'fontsize',10,'box','on','linewidth',1);
    errorbar(Y1,Y2,Y2Err,Y2Err,Y1Err,Y1Err,'o','MarkerFaceColor',co(1,:),...
        'Color',co(1,:)*.5,'linewidth',1,'markersize',8);
    xlabel('xdt1 y-position (mV/V)')
    ylabel('xdt2 y-position (mV/V)')
   
