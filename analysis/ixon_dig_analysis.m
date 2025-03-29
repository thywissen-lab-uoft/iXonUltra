% ixon_main_dig.m
%
% Author : CF Fujiwara
%
% This script is the primary analysis code for ixon images which are single
% plane.

% Display this filename
disp(repmat('-',2,60));disp(mfilename);disp(repmat('-',2,60)); 
% Add all subdirectories for this m file
curpath = fileparts(mfilename('fullpath'));
addpath(curpath);addpath(genpath(curpath));
%% Close all non GUI figures
% Close all figures without the GUI tag.
figs=get(groot,'Children');
disp('Closing all non GUI figures.');
for kk=1:length(figs)
   if ~isequal(figs(kk).Tag,'GUI')
        disp(['Closing figure ' num2str(figs(kk).Number) ' ' figs(kk).Name]);
        close(figs(kk)) 
   end
end
disp(' ');
%% Select image directory
    
if ~exist('dig_auto_file');dig_auto_file = 1;end

if dig_auto_file
    dialog_title='Select digdata';       
    [filename,dig_imgdir,b]=uigetfile(fullfile(ixon_getDayDir,'*.mat'),dialog_title);
    filename = fullfile(dig_imgdir,filename);    
    if  b == 0
        disp('Canceling.');    
        return; 
    end
end
%% Load the Data
clear digdata
tic
fprintf('Loading digdata ...');
digdata = load(filename);
disp([' done (' num2str(toc,'%.2f') 's)']);

%% Initialize Options

dig_opts = struct;
dig_opts.Quality = 'auto';
dig_opts.saveDir=dig_imgdir;
dig_opts.FigLabel=digdata.SourceDirectory{1};

%% Analysis Variable
% This section of code chooses the variable to plot against for aggregate
% plots.  The chosen variable MUST match a variable provided in the params
% field of the .mat file. The unit has no tangibile affect and only affects
% display properties.

% Choose what kind of variable to plot against (sequencer/camera)
dig_opts.varType        = 'param';          % always select 'param' for now 
dig_opts.autoXVar       = 1;                % Auto detect changing variable?
dig_opts.autoUnit       = 1;                % Auto detect unit for variable?
dig_opts.xVar           = 'conductivity_mod_time';  % Variable Name
dig_opts.overrideUnit   = 'V';              % If ixon_autoUnit=0, use this
dig_opts.doSave         = 1;                % Save Analysis?
autoVar_Ignore = {'f_offset','piezo_offset'};

%% Flags

% Recenter all binned data to have same limits
dig_doShowCloud                         = 1;
dig_doShowCloudAnimate                  = 1;
dig_standardAnalysis                    = 1;
dig_ac_conductivity_fit                 = 0;
dig_quench_conductivity_fit             = 0;
dig_doRadialAnalysis                        = 0; % has issues,obsolete
dig_doRadialSkewAnalysis                    = 0; % has issues,obsolete

dig_doRadialAnalysis2                   = 1;
dig_doFidelity                          = 1;

do_qpd_analysis                         = 1;

%% QPD Analysis
if do_qpd_analysis
    try
    P=[digdata.Params];
    D=[P.ExecutionDate];
    dig_opts.X = [digdata.X];
    dig_opts.xVar = digdata.xVar;
    
    dig_opts.Frequency = unique([digdata.Params.conductivity_mod_freq]);
    dig_opts.RampTime  = unique([digdata.Params.conductivity_mod_ramp_time]);

   [figs,output] = qpd_main(D,dig_opts) ;
      if dig_opts.doSave
          for kk=1:length(figs)
              if ~isempty(figs{kk})
                ixon_saveFigure2(figs{kk},...
                 figs{kk}.Name,dig_opts);  
              end
          end
           try if ~exist(dig_opts.saveDir,'dir');mkdir(dig_opts.saveDir);end;end
            filename = fullfile(dig_opts.saveDir,'qpd.mat');
            disp(['Saving ' filename ' ...']);
            save(filename, '-struct','output');      
      end  
    catch ME
       warning(getReport(ME),'extended','hyperlinks','on');
       disp('DID YOU MAKE SURE TO ADD THE AUXLIARY_ANALYSIS GIT REPOSITORY IS LOADE?');
    end
end

%% Basic Analysis
digdata = dig_basic(digdata);

%% Show Cloud

if dig_doShowCloud
    opts = dig_opts;
    opts.doAnimate = 1;
    opts.ROI = 'max';
    opts.filename = 'dig_animateFull';
    hF_digCloud = dig_showCloud(digdata,opts);
    if dig_opts.doSave
        ixon_saveFigure2(hF_digCloud,...
         'dig_average',dig_opts);  
    end
end

%% Standard Analysis

if dig_standardAnalysis
    hF_digStandard = dig_showStandardAnalysis(digdata,opts);
    if dig_opts.doSave
        ixon_saveFigure2(hF_digStandard,...
         'dig_standard',dig_opts);  
    end
end

%% Fidelity
if dig_doFidelity && size(digdata.Zdig,4)==2
    % CF Needs to finish writing this
    digdata = dig_Fidelity(digdata);
    hF_FidelityMap = dig_showFidelityMap(digdata,dig_opts);
    if dig_opts.doSave
       
        if ~isempty(hF_FidelityMap)
                ixon_saveFigure2(hF_FidelityMap,...
                    hF_FidelityMap.Name,dig_opts);  
        end           
    end
end

%% Radial Analysis Better

if dig_doRadialAnalysis2
    [hF,dig_radial_data] = dig_radialAnalysis(digdata);
    if dig_opts.doSave
        ixon_saveFigure2(hF,'dig_radial',dig_opts);  
    end

    try if ~exist(dig_opts.saveDir,'dir');mkdir(dig_opts.saveDir);end;end
    filename = fullfile(dig_opts.saveDir,'dig_radial_data.mat');
    disp(['Saving ' filename ' ...']);
    save(filename, '-struct','dig_radial_data'); 
end

%% Radial Analysis
% CJF Thinks is becoming obsolte

if dig_doRadialAnalysis
    opts = dig_opts;   
    opts.BinStep = 3;                   % delta r bin size    
    opts.useAverageCenter = 0;

    if isequal(digdata.xVar,'conductivity_mod_time')
        opts.useAverageCenter = 1;
    end

    digdata = dig_compute_radial(digdata,opts);     % Compute radial profile for all images


    
    opts.rMaxShow = 80;                 % max r to plot
    opts.nMaxShow = 1;               % max density to plot
    opts.showDevParametrization  = 0;   % show standard deviation?
    
    % Show radial profiles
    hFs=dig_showRadialProfile(digdata,opts);     
     if dig_opts.doSave
         for kk=1:length(hFs)
             ixon_saveFigure2(hFs(kk),...
                ['dig_radial_profile_' num2str(kk)],dig_opts);  
         end
     end        
     
      opts.ForceAverage = 0;
     if isequal(digdata.xVar,'ExecutionDate')
        opts.ForceAverage = 1; 
     end
    if isequal(digdata.xVar,'conductivity_mod_time')
        opts.ForceAverage = 1;
    end

%     opts.ForceAverage = 0;
    opts.doAnimate = 1;
    
    [hFs_radial,hF2_radial] = dig_radialAnalysis_average_images(digdata,opts);
    if dig_opts.doSave
       for kk=1:length(hFs_radial)
           ixon_saveFigure2(hFs_radial(kk),...
                ['dig_radial_' num2str(kk)],dig_opts);   
       end
       
       ixon_saveFigure2(hF2_radial,...
                ['dig_filling'],dig_opts);   
       
    end
    
    doExportforDrut=0;
    if doExportforDrut    
       % Constants
        h = 6.62607015*10^-34;
        m = 39.96399848*1.66053906660*10^-27;
        lam = (1054*10^-9*1054*10^-9*1064*10^-9)^(1/3);
        aL = lam/2;

        % 120mW XDT + 2.5ER request lattice trap frequencies - calibrated
        % 02/29/24 (reanalyzed 04/11/24)
        omega_x = 2*pi*67.3;
        omega_y = 2*pi*60.2;
        omega_bar = (omega_x*omega_y).^(1/2);

        % Harmonic potential in Hz
        PotentialVector = 0.5*m*omega_bar^2*aL^2.*rVec.^2/h;

        kappa = 0.5*m*omega_bar^2*aL^2/h;

        % Assign to output
        out = struct;
        out.Images                         = digdata.Z;
        out.SiteVector1                    = digdata.n1;
        out.SiteVector2                    = digdata.n2;            
        out.RadialVector                   = digdata.r;
        out.RadialOccupation               = digdata.Zr;
        out.RadialPotential                = digdata.r.^2*kappa;                  

        try if ~exist(dig_opts.saveDir,'dir');mkdir(dig_opts.saveDir);end;end
        filename = fullfile(dig_opts.saveDir,'dig_radial_data.mat');
        disp(['Saving ' filename ' ...']);
        save(filename, '-struct','dig_radial_data');
    end    

end

%% Radial Skew Analysis

if  dig_doRadialSkewAnalysis && isequal(digdata.xVar,'ExecutionDate')    
    [digdata] = dig_compute_radial_skew(digdata,opts);
    [hFs_radial_skew] = dig_radialAnalysis_average_images_skew(digdata,opts);
    if dig_opts.doSave
       for kk=1:length(hFs_radial_skew)
           ixon_saveFigure2(hFs_radial_skew(kk),...
                ['dig_radial_skew_' num2str(kk)],dig_opts);   
       end
    end
end

%% Conductivity Analysis

if dig_ac_conductivity_fit
opts = dig_opts;

opts.QPD_phi = mean([output.QPD_Modulation.Phi1 output.QPD_Modulation.Phi2]);
% opts.QPD_phi = 0;
[hF_conductivity,conductivity_data] = dig_ac_conductivity(digdata,opts);
    if dig_opts.doSave
        ixon_saveFigure2(hF_conductivity,...
         'dig_conductivity',dig_opts);  
    end
    try if ~exist(dig_opts.saveDir,'dir');mkdir(dig_opts.saveDir);end;end
    filename = fullfile(dig_opts.saveDir,'conductivity_data.mat');
    disp(['Saving ' filename ' ...']);
    save(filename, '-struct','conductivity_data');
% keyboard

end

if dig_quench_conductivity_fit
opts = dig_opts;

[hF_quench,quench_data] = dig_quench_conductivity(digdata,opts);
    if dig_opts.doSave
        ixon_saveFigure2(hF_quench,...
         'dig_quench',dig_opts);  
    end
    try if ~exist(dig_opts.saveDir,'dir');mkdir(dig_opts.saveDir);end;end
    filename = fullfile(dig_opts.saveDir,'quench_data.mat');
    disp(['Saving ' filename ' ...']);
    save(filename, '-struct','quench_data');


end

