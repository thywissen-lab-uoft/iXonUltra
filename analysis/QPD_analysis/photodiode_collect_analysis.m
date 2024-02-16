function [pd_summary,pd_data] = photodiode_collect_analysis(qpdfiles,opts)

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
    pd_data(nn) = photodiode_analyze(qpdfiles(nn));
end

%% Assign analyzed data to output

pd_out.qpd_data = pd_data;

%% Compute average parameters for all runs

% Take average and standard deviation of fit values during first 50ms
    % of modulation

pd_out.mod_params.X1.A     = mean([pd_data.modfit_X1_A]);
pd_out.mod_params.X1.A_err = std([pd_data.modfit_X1_A]);
pd_out.mod_params.X1.B     = mean([pd_data.modfit_X1_B]);
pd_out.mod_params.X1.B_err = std([pd_data.modfit_X1_B]);
pd_out.mod_params.X1.C     = mean([pd_data.modfit_X1_C]);
pd_out.mod_params.X1.C_err = std([pd_data.modfit_X1_C]);
pd_out.mod_params.X1.T     = mean([pd_data.modfit_X1_T]);
pd_out.mod_params.X1.T_err = std([pd_data.modfit_X1_T]);

pd_out.mod_params.X2.A     = mean([pd_data.modfit_X2_A]);
pd_out.mod_params.X2.A_err = std([pd_data.modfit_X2_A]);
pd_out.mod_params.X2.B     = mean([pd_data.modfit_X2_B]);
pd_out.mod_params.X2.B_err = std([pd_data.modfit_X2_B]);
pd_out.mod_params.X2.C     = mean([pd_data.modfit_X2_C]);
pd_out.mod_params.X2.C_err = std([pd_data.modfit_X2_C]);
pd_out.mod_params.X2.T     = mean([pd_data.modfit_X2_T]);
pd_out.mod_params.X2.T_err = std([pd_data.modfit_X2_T]);

% Take average and standard deviation of average powers during -350ms to 50ms
pd_out.powers.ODT1      = mean([pd_data.ODT1_ave]);
pd_out.powers.ODT1_err  = std([pd_data.ODT1_ave]);
pd_out.powers.ODT2      = mean([pd_data.ODT2_ave]);
pd_out.powers.ODT2_err  = std([pd_data.ODT2_ave]);
pd_out.powers.XLatt     = mean([pd_data.XLatt_ave]);
pd_out.powers.XLatt_err = std([pd_data.XLatt_ave]);
pd_out.powers.YLatt     = mean([pd_data.YLatt_ave]);
pd_out.powers.YLatt_err = std([pd_data.YLatt_ave]);
pd_out.powers.ZLatt     = mean([pd_data.ZLatt_ave]);
pd_out.powers.ZLatt_err = std([pd_data.ZLatt_ave]);


end