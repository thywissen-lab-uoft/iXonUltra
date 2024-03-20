function output = photodiode_analyze(qpd_filenames,params)

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analyzes a single QPD file
% Outputs average powers, fit parameters, raw channels with time shifted to
% start of modulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear qpd_out;

disp(['Analysing ' num2str(length(qpd_filenames)) ' qpd files...']);

for nn=1:length(qpd_filenames)
    fprintf(['(' num2str(nn) '/' num2str(length(qpd_filenames)) ') ']);

    %% Initialize data struct and load QPD data

    qpd_single = load(qpd_filenames{nn});
    
    [~,b,~]=fileparts(qpd_filenames{nn});
    fprintf([b ' ']);

    %% Shift the time such that modulation start at t=0ms

    % Index of modulation start and convert to ms
    mod_start = max(find(diff(qpd_single.data(:,10))==1));

    % Shift time by modulation start time
    t = (qpd_single.t - qpd_single.t(mod_start))*1e3;

    %% Save all traces with timeshift

    output(nn).t       = t;
    output(nn).X1      = qpd_single.data(:,1);
    output(nn).Y1      = qpd_single.data(:,2);
    output(nn).SUM1    = qpd_single.data(:,3);
    output(nn).X2      = qpd_single.data(:,4);
    output(nn).Y2      = qpd_single.data(:,5);
    output(nn).SUM2    = qpd_single.data(:,6);
    output(nn).XLATT   = qpd_single.data(:,7);
    output(nn).YLATT   = qpd_single.data(:,8);
    output(nn).ZLATT   = qpd_single.data(:,9);
    output(nn).TRIG    = qpd_single.data(:,10);
    
    %% Analyze the lattice load traces
    fprintf('lattice...');
    
    % Period of 60Hz in ms
    T60 = 1e3/60;
    
    % Indeces before modulation
    ii = t<=0;

    % Get the lattice photodiode information
    xlatt = output(nn).XLATT;
    ylatt = output(nn).YLATT;
    zlatt = output(nn).ZLATT;

    % Get the photodiode information before loading
    tpre  = t(ii);
    xlatt = xlatt(ii);
    ylatt = ylatt(ii);
    zlatt = zlatt(ii);
    
    % How many 60 Hz periods to measure the zero level
    tzero = T60*6;    
    
    iL = find(tpre<=tzero);
    
    % Measure the zero level mean
    xlatt_L_mean = mean(xlatt(iL));
    ylatt_L_mean = mean(ylatt(iL));
    zlatt_L_mean = mean(zlatt(iL));

    % Measure the zero level standard deviation
    xlatt_L_std = std(xlatt(iL));
    ylatt_L_std = std(ylatt(iL));
    zlatt_L_std = std(zlatt(iL));
    
    % How many 60 Hz periods to measure the high level
    thigh = T60*6;    
    
    iH = find([tpre<=0].*[tpre>=-thigh]);
    
    % Measure the zero level mean
    xlatt_H_mean = mean(xlatt(iH));
    ylatt_H_mean = mean(ylatt(iH));
    zlatt_H_mean = mean(zlatt(iH));

    % Measure the zero level standard deviation
    xlatt_H_std = std(xlatt(iH));
    ylatt_H_std = std(ylatt(iH));
    zlatt_H_std = std(zlatt(iH));
    
    % Measure the x lattice delta
    xlatt_val = xlatt_H_mean-xlatt_L_mean;
    xlatt_err = sqrt(xlatt_H_std.^2+xlatt_L_std.^2);

    % Measure the y lattice delta
    ylatt_val = ylatt_H_mean-ylatt_L_mean;
    ylatt_err = sqrt(ylatt_H_std.^2+ylatt_L_std.^2);
    
    % Measure the z lattice delta
    zlatt_val = zlatt_H_mean-zlatt_L_mean;
    zlatt_err = sqrt(zlatt_H_std.^2+zlatt_L_std.^2);
    
    output(nn).xlatt_val = xlatt_val;
    output(nn).xlatt_err = xlatt_err;
    output(nn).ylatt_val = ylatt_val;
    output(nn).ylatt_err = ylatt_err;
    output(nn).zlatt_val = zlatt_val;
    output(nn).zlatt_err = zlatt_err;
    %% Average the power over the the 200ms before the ramp --> first 50ms of mod.
   fprintf('odt powers...');
   
    % Select time to average over - right now choosing 200 ms before ramp
    % time, and first 50 ms of modulation
    % modulation 

    start_ave = find(round(t,1)==-350);
    end_ave   = find(round(t,1)==50);

    output(nn).t_ave = t(start_ave:end_ave);

    % Average the powers over the averaging window

    output(nn).ODT1_ave = mean(qpd_single.data(start_ave:end_ave,3));
    output(nn).ODT1_std = std(qpd_single.data(start_ave:end_ave,3));

    output(nn).ODT2_ave = mean(qpd_single.data(start_ave:end_ave,6));
    output(nn).ODT2_std = std(qpd_single.data(start_ave:end_ave,6));

    output(nn).XLatt_ave = mean(qpd_single.data(start_ave:end_ave,7));
    output(nn).XLatt_std = std(qpd_single.data(start_ave:end_ave,7));

    output(nn).YLatt_ave = mean(qpd_single.data(start_ave:end_ave,8));
    output(nn).YLatt_std = std(qpd_single.data(start_ave:end_ave,8));

    output(nn).ZLatt_ave = mean(qpd_single.data(start_ave:end_ave,9));
    output(nn).ZLatt_std = std(qpd_single.data(start_ave:end_ave,9));

    %% Fit the first 50 ms of modulation
    fprintf('odt modulation...');

    % Select the first 50 ms of modulation
    mod50 = find(round(t,1)==50);
    t_fitmod = t(mod_start:mod50);

    % Normalize the data by the sum
    X1_fitmod = qpd_single.data(mod_start:mod50,1)./qpd_single.data(mod_start:mod50,3);
    X2_fitmod = qpd_single.data(mod_start:mod50,4)./qpd_single.data(mod_start:mod50,6);

    % Basic sin func to fit
    myfunc = @(A,B,C,T,t) A*sin(2*pi*t/T + B) + C;    

    myfit = fittype(@(A,B,C,T,t) myfunc(A,B,C,T,t),'independent',{'t'},...
        'coefficients',{'A','B','C','T'});

    Agp = (max(X1_fitmod)-min(X1_fitmod))/2;
    Cgp = mean(X1_fitmod);

%     [pks,plocs] = findpeaks(smoothdata(X1_fitmod,"gaussian",10),'MinPeakProminence',Agp/4);
%     [trs,tlocs] = findpeaks(-smoothdata(X1_fitmod,"gaussian",10),'MinPeakProminence',Agp/4);
%     npeaks = min(length(plocs),length(tlocs));
%     Tgp = 2*abs(mean(t_fitmod(tlocs(1:npeaks))-t_fitmod(plocs(1:npeaks))));
%     Tgp = 1/0.06;
%     Tgp = 1/(0.145);

    if isfield(params(nn),'conductivity_mod_freq')
        Tgp = 1e3/(params(nn).conductivity_mod_freq);
    else
        Tgp = 1/0.06;
    end
    
    Bgp =  mod(-2*pi*(50)/Tgp,2*pi);
    opt = fitoptions(myfit);
    opt.StartPoint = [Agp Bgp Cgp Tgp];
    opt.Lower  = [0 -10*pi Cgp-20 0];
    fout_X1 = fit(t_fitmod',X1_fitmod,myfit,opt);  

    Agp = (max(X2_fitmod)-min(X2_fitmod))/2;
    Cgp = mean(X2_fitmod);
    Bgp =  mod(2*pi*(50)/Tgp,2*pi);
    opt = fitoptions(myfit);
    opt.StartPoint = [Agp Bgp Cgp Tgp];
    opt.Lower  = [0 -10*pi Cgp-20 0];
    fout_X2 = fit(t_fitmod',X2_fitmod,myfit,opt);
%%

    output(nn).modfit_t      = t_fitmod;
    output(nn).modfit_normX1data = X1_fitmod;
    output(nn).modfit_normX2data = X2_fitmod;

    output(nn).modfit_X1     = fout_X1;
    output(nn).modfit_X2     = fout_X2;

    output(nn).modfit_X1_A   = fout_X1.A;
    output(nn).modfit_X2_A   = fout_X2.A;
    output(nn).modfit_X1_B   = mod(fout_X1.B,2*pi);
    output(nn).modfit_X2_B   = mod(fout_X2.B,2*pi);
    output(nn).modfit_X1_C   = fout_X1.C;
    output(nn).modfit_X2_C   = fout_X2.C;
    output(nn).modfit_X1_T   = fout_X1.T;
    output(nn).modfit_X2_T   = fout_X2.T;

disp('done');
end

end
