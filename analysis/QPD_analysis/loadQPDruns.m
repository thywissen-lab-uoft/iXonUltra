function [ixondata,qpd_out] = loadQPDruns(ixondata,qpdfiles,opts)

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analyzes list of QPD filenames 
% Outputs average powers, fit parameters, raw channels with time shifted to
% start of modulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize options

if nargin~=2
    opts.doPlot = 1;
    opts.doSave = 0;
end

% - should plotting be separated from analysis?

%% Analyze the list of QPD files

% Initialize data struct
% qpd_out = struct;

for nn=1:length(qpdfiles)

    % Analyze single QPD file
    qpd_data(nn) = loadQPD(qpdfiles(nn));
    ixondata(nn).qpd_data = qpd_data(nn);
end

%% Assign analyzed data to output

qpd_out.qpd_data = qpd_data;

%% Compute average parameters for all runs

% Take average and standard deviation of fit values during first 50ms
    % of modulation

qpd_out.mod_params.X1.A     = mean([qpd_data.modfit_X1_A]);
qpd_out.mod_params.X1.A_err = std([qpd_data.modfit_X1_A]);
qpd_out.mod_params.X1.B     = mean([qpd_data.modfit_X1_B]);
qpd_out.mod_params.X1.B_err = std([qpd_data.modfit_X1_B]);
qpd_out.mod_params.X1.C     = mean([qpd_data.modfit_X1_C]);
qpd_out.mod_params.X1.C_err = std([qpd_data.modfit_X1_C]);
qpd_out.mod_params.X1.T     = mean([qpd_data.modfit_X1_T]);
qpd_out.mod_params.X1.T_err = std([qpd_data.modfit_X1_T]);

qpd_out.mod_params.X2.A     = mean([qpd_data.modfit_X2_A]);
qpd_out.mod_params.X2.A_err = std([qpd_data.modfit_X2_A]);
qpd_out.mod_params.X2.B     = mean([qpd_data.modfit_X2_B]);
qpd_out.mod_params.X2.B_err = std([qpd_data.modfit_X2_B]);
qpd_out.mod_params.X2.C     = mean([qpd_data.modfit_X2_C]);
qpd_out.mod_params.X2.C_err = std([qpd_data.modfit_X2_C]);
qpd_out.mod_params.X2.T     = mean([qpd_data.modfit_X2_T]);
qpd_out.mod_params.X2.T_err = std([qpd_data.modfit_X2_T]);


% Take average and standard deviation of average powers during -350ms to 50ms

qpd_out.powers.ODT1      = mean([qpd_data.ODT1_ave]);
qpd_out.powers.ODT1_err  = std([qpd_data.ODT1_ave]);
qpd_out.powers.ODT2      = mean([qpd_data.ODT2_ave]);
qpd_out.powers.ODT2_err  = std([qpd_data.ODT2_ave]);
qpd_out.powers.XLatt     = mean([qpd_data.XLatt_ave]);
qpd_out.powers.XLatt_err = std([qpd_data.XLatt_ave]);
qpd_out.powers.YLatt     = mean([qpd_data.YLatt_ave]);
qpd_out.powers.YLatt_err = std([qpd_data.YLatt_ave]);
qpd_out.powers.ZLatt     = mean([qpd_data.ZLatt_ave]);
qpd_out.powers.ZLatt_err = std([qpd_data.ZLatt_ave]);


end