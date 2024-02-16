function qpd_out = photodiode_analyze(qpd_filename)

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analyzes a single QPD file
% Outputs average powers, fit parameters, raw channels with time shifted to
% start of modulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for nn=1:length(qpd_filename)
    %% Initialize data struct and load QPD data

    qpd_single = load(qpd_filename);

    %% Shift the time such that modulation start at t=0ms

    % Index of modulation start and convert to ms
    mod_start = max(find(diff(qpd_single.data(:,10))==1));

    % Shift time by modulation start time
    t = (qpd_single.t - qpd_single.t(mod_start))*1e3;

    %% Save all of the traces with the correct shifted time

    qpd_out(nn).t       = t;
    qpd_out(nn).X1      = qpd_single.data(:,1);
    qpd_out(nn).Y1      = qpd_single.data(:,2);
    qpd_out(nn).SUM1    = qpd_single.data(:,3);
    qpd_out(nn).X2      = qpd_single.data(:,4);
    qpd_out(nn).Y2      = qpd_single.data(:,5);
    qpd_out(nn).SUM2    = qpd_single.data(:,6);
    qpd_out(nn).XLATT   = qpd_single.data(:,7);
    qpd_out(nn).YLATT   = qpd_single.data(:,8);
    qpd_out(nn).ZLATT   = qpd_single.data(:,9);
    qpd_out(nn).TRIG    = qpd_single.data(:,10);

    %% Average the power over the the 200ms before the ramp --> first 50ms of mod.

    % Select time to average over - right now choosing 200 ms before ramp
    % time, and first 50 ms of modulation
    % modulation 

    start_ave = find(round(t,1)==-350);
    end_ave   = find(round(t,1)==50);

    qpd_out(nn).t_ave = t(start_ave:end_ave);

    % Average the powers over the averaging window

    qpd_out(nn).ODT1_ave = mean(qpd_single.data(start_ave:end_ave,3));
    qpd_out(nn).ODT1_std = std(qpd_single.data(start_ave:end_ave,3));

    qpd_out(nn).ODT2_ave = mean(qpd_single.data(start_ave:end_ave,6));
    qpd_out(nn).ODT2_std = std(qpd_single.data(start_ave:end_ave,6));

    qpd_out(nn).XLatt_ave = mean(qpd_single.data(start_ave:end_ave,7));
    qpd_out(nn).XLatt_std = std(qpd_single.data(start_ave:end_ave,7));

    qpd_out(nn).YLatt_ave = mean(qpd_single.data(start_ave:end_ave,8));
    qpd_out(nn).YLatt_std = std(qpd_single.data(start_ave:end_ave,8));

    qpd_out(nn).ZLatt_ave = mean(qpd_single.data(start_ave:end_ave,9));
    qpd_out(nn).ZLatt_std = std(qpd_single.data(start_ave:end_ave,9));

    %% Fit the first 50 ms of modulation

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
    Tgp = 1/(0.06);
    Bgp =  mod(2*pi*(150)/Tgp,2*pi);
    opt = fitoptions(myfit);
    opt.StartPoint = [Agp Bgp Cgp Tgp];
    opt.Lower  = [0 -10*pi Cgp-20 0];
    fout_X1 = fit(t_fitmod',X1_fitmod,myfit,opt);

    Agp = (max(X2_fitmod)-min(X2_fitmod))/2;
    Cgp = mean(X2_fitmod);
    Bgp =  mod(2*pi*(150)/Tgp,2*pi);
    opt = fitoptions(myfit);
    opt.StartPoint = [Agp Bgp Cgp Tgp];
    opt.Lower  = [0 -10*pi Cgp-20 0];
    fout_X2 = fit(t_fitmod',X2_fitmod,myfit,opt);

    qpd_out(nn).modfit_t      = t_fitmod;
    qpd_out(nn).modfit_normX1data = X1_fitmod;
    qpd_out(nn).modfit_normX2data = X2_fitmod;

    qpd_out(nn).modfit_X1     = fout_X1;
    qpd_out(nn).modfit_X2     = fout_X2;

    qpd_out(nn).modfit_X1_A   = fout_X1.A;
    qpd_out(nn).modfit_X2_A   = fout_X2.A;
    qpd_out(nn).modfit_X1_B   = mod(fout_X1.B,2*pi);
    qpd_out(nn).modfit_X2_B   = mod(fout_X2.B,2*pi);
    qpd_out(nn).modfit_X1_C   = fout_X1.C;
    qpd_out(nn).modfit_X2_C   = fout_X2.C;
    qpd_out(nn).modfit_X1_T   = fout_X1.T;
    qpd_out(nn).modfit_X2_T   = fout_X2.T;
end

end