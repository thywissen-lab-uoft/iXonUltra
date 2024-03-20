% pd_analysis.m
%
% Author : CF Fujiwara
%

% Display this filename
disp(repmat('-',1,60));disp(repmat('-',1,60));    
disp(['Calling ' mfilename '.m']);
disp(repmat('-',1,60));disp(repmat('-',1,60));    


%% Flags

PD_do_lattice_load = 1;
PD_do_save = 1;
PD_do_lattice_fluor = 1;
PD_do_odt_power = 1;
PD_do_odt_modulation = 1;

%% Initialize Data
for nn=1:length(pd_data)
    
    % Grab data and convert to mV and ms
    data = pd_data(nn).Data.data;
    data(:,1:9) = data(:,1:9)*1e3;    
    t = 1e3*pd_data(nn).Data.t;  
    % Index of modulation start
    i0 = max(find(diff(data(:,10))==1));
    % Shift time by modulation start time
    t = t - t(i0);
    pd_data(nn).t = t;
    pd_data(nn).v = data;
end


%% Lattice Load Analysis
if PD_do_lattice_load
    opts = struct;
    opts.FigLabel = pd_opts.FigLabel;
    opts.XCalibration = 219.74/(950-14); % Er/mV
    opts.XCalibrationUncertainty = 0.563/(950-14); % Er/mV
    opts.XCalibrationString = '2023-03-11'; % date please

    opts.YCalibration = 196.15/(1510-14.4); % Er/mV
    opts.YCalibrationUncertainty = 1.15/(1510-14.4); % Er/mV
    opts.YCalibrationString = '2023-03-11'; % date please
    
    opts.ZCalibration = 205.89/(1020-13.2); % Er/mV
    opts.ZCalibrationUncertainty = 9.13/(1020-13.2); % Er/mV
    opts.ZCalibrationString = '2023-03-11'; % date please
    
    [pd_data,lattice_summary,hF_lattice_load] = pd_lattice_load(pd_data,opts);
    
    if PD_do_save
        ixon_saveFigure2(hF_lattice_load,'pd_lattice_load',pd_opts);
    end
end

%% Lattice Fluoresence Analysis
if PD_do_lattice_fluor
    opts = struct;
    opts.FigLabel = pd_opts.FigLabel;    
    hF_lattice_fluor = pd_lattice_fluor(pd_data,opts);
    
    if PD_do_save
        ixon_saveFigure2(hF_lattice_fluor,'pd_lattice_fluor',pd_opts);
    end
end

%% ODT Power Analysis
if PD_do_odt_power
    opts = struct;
    opts.FigLabel = pd_opts.FigLabel;    
    [pd_data,odt_summary,hF_odt_power] = pd_odt_power(pd_data,opts);
    
    if PD_do_save
        ixon_saveFigure2(hF_odt_power,'pd_odt_power',pd_opts);
    end
end

%% ODT Modulation Analysis
if PD_do_odt_modulation
    opts = struct;
    opts.FigLabel = pd_opts.FigLabel;    
    
    opts.ODT1Calibration  = [94 162.5]; % um/V/V xlattice, um/V/V ylattice
    opts.ODT2Calibration  = [995 -617]; % um/V/V xlattice, um/V/V ylattice
    opts.ODTCalibrationString = '2024-02-13';

    [pd_data,modulation_summary,hF_odt_modulation] = pd_odt_modulation(pd_data,opts);
    
    if PD_do_save
        ixon_saveFigure2(hF_odt_modulation,'pd_odt_modulation',pd_opts);
    end    
end