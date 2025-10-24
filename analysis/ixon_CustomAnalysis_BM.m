%% ixon_CustomAnalysis_BM.m
%
% This script runs customized analysis on the box counts or gaussian fits.
% In particular, it enables you to customized how the analysis is
% performed.

%% Custom Analysis Data Source

% Select the data source

data_source = 'box';
% data_source = 'gauss';

switch data_source
    case 'box'
        src_data = ixon_boxdata;
    case 'gauss'
        src_data = ixon_gaussdata;
end

%Select which image number you want to analyze
imNum = 1;

%% Y Data Select Flags

y_Lbl = {};

%%%%%% Total Number
y_Lbl{01}       = 'N_T';
%%%%%% Relative Ground Band
y_Lbl{02}       = 'N_g/N_T';
%%%%%% Relative Px Band
y_Lbl{03}       = 'N_{P_x}/N_T';
%%%%%% Relative Py Band
y_Lbl{04}       = 'N_{P_y}/N_T';

p_inds = [1,2,3,4];

%% Fit Flags (Fitting uses customFit.m, same as 

FitFlags = struct;

% Exponential Decay
FitFlags.expdecay = 0;               % Exponential Decay no offset
FitFlags.expdecayoffset = 0;        % Exponential Decay w/offset


FitFlags.T2exp=0;


FitFlags.Rabi_oscillation = 0;
FitFlags.Rabi_oscillation2 = 0;
FitFlags.NGaussPeak=0;

FitFlags.gauss_single=0;
FitFlags.gauss_4=0;
FitFlags.gauss_neg_double=0;
FitFlags.gauss_neg_single=0;
FitFlags.gauss_double = 0;
FitFlags.gauss_triple = 0;
 
FitFlags.lorentz_neg_single=0;    
FitFlags.lorentz_neg_double=0;  

FitFlags.lorentz_single=0;
FitFlags.lorentz_double=0;    
FitFlags.lorentz_triple=0;    

FitFlags.lorentz_asym_single= 0;
FitFlags.lorentz_asym_double= 0;

FitFlags.fit_lorentz_assymetric_4=0;

%% X Data

doCustomX = 0; %Custom X needs to be fixed for ixon analysis

custom_data_bm = struct;
custom_data_bm.Source = src_data;

% Assign atom number
custom_data_bm.Natoms(:,:,1) = src_data.N(:,imNum,:); 

% Get the default X data;
custom_data_bm.X = src_data.X;
custom_data_bm.XVar = src_data.xVar;
custom_data_bm.XStr = src_data.xVar;
custom_data_bm.XUnit = src_data.Units(1).(src_data.xVar);

if doCustomX
    % Select mF states
    mF1 = -7/2;
    mF2 = -9/2;
     
    % mF1 = -7/2;
    % mF2 = -5/2;

   Bfb   = src_data.Params(1).HF_FeshValue_Initial_Lattice;
%     Bfb   = src_data.Params(1).HF_FeshValue_Spectroscopy;
%     Bfb   = src_data.Params(1).HF_FeshValue_Final_Lattice;
    Bshim = src_data.Params(1).HF_zshim_Initial_Lattice*2.35;
%     Boff  = 0.11;
    Boff  = 0.11;
    B = Bfb + Bshim + Boff;
%     B = 199.5 + 0 + 0.11; 

    
    % Transition Energy
    x0 = abs((BreitRabiK(B,9/2,mF1)-BreitRabiK(B,9/2,mF2)))/6.6260755e-34/1E6; 

    % Convert x variable into transition energy
    switch src_data.xVar
        case 'Raman_AOM3_freq'
            X=custom_data_bm.X;
            X = 2*X - 80;
            X = X - x0;   
            X = X*1e3;
            xstr=['frequency - ' num2str(round(abs(x0),4))  ' MHz (kHz)']; 
            xunit = 'kHz';
        case 'Pulse_Time'
            X=custom_data_bm.X;
            xstr='pulse time (ms)';    
            xunit = 'ms';
       case 'rf_rabi_time_HF'
            X=custom_data_bm.X;
            xstr='pulse time (ms)';    
            xunit = 'ms';
        case 'rf_freq_HF'
            X=custom_data_bm.X;
            X = X - x0;   
            X = X*1e3;
            xstr=['frequency - ' num2str(round(abs(x0),4))  ' MHz (kHz)']; 
            xunit = 'kHz';
        case 'rf_rabi_freq_HF'
            X=custom_data_bm.X;
            X = X - x0;   
            X = X*1e3;
            xstr=['frequency - ' num2str(round(abs(x0),4))  ' MHz (kHz)']; 
            xunit = 'kHz';
        case 'rf_tof_freq'
          X=custom_data_bm.X;
          B = src_data.Params(1).HF_FeshValue_Final_Lattice + 0 + 0.11; 

    
            % Transition Energy
            x0 = abs((BreitRabiK(B,9/2,mF1)-BreitRabiK(B,9/2,mF2)))/6.6260755e-34/1E6; 
            X = X - x0;   
            X = X*1e3;
            xstr=['frequency - ' num2str(round(abs(x0),4))  ' MHz (kHz)'];  
            xunit = 'kHz';
        otherwise
            X = custom_data_bm.X;
            xstr = pco_xVar;
            xunit = 'ms';
    end 
    
    % Assign outputs
    custom_data_bm.x0 = x0;
    custom_data_bm.X = X;
    custom_data_bm.XStr = xstr;        
    custom_data_bm.XLabel = xstr;        
    custom_data_bm.XUnit = xunit;
end    

%% Y Data

NT = sum(custom_data_bm.Natoms(:,2:6),2);

Ng = custom_data_bm.Natoms(:,2);
Npy = custom_data_bm.Natoms(:,3) + custom_data_bm.Natoms(:,4);
Npx = custom_data_bm.Natoms(:,5) + custom_data_bm.Natoms(:,6);

Y = struct;

%%%%%%% Total Number
Y(1).YName      = y_Lbl{01};
Y(1).FigName    = 'ixon_BM_custom_NTot';
Y(1).Y          = NT;

Y(2).YName      = y_Lbl{02};
Y(2).FigName    = 'ixon_BM_custom_Ng_rel';
Y(2).Y          = Ng./NT;

Y(3).YName      = y_Lbl{03};
Y(3).FigName    = 'ixon_BM_custom_NPx_rel';
Y(3).Y          = Npx./NT;

Y(4).YName      = y_Lbl{04};
Y(4).FigName    = 'ixon_BM_custom_NPy_rel';
Y(4).Y          = Npy./NT;


% Assign to output
custom_data_bm.Y = Y;
custom_data_bm.YLabel = {Y.YName};        

%% Plot it 
bm_custom_opts=struct;
fouts={};
names={};
clear hFs

% X = X - 391016.821;
X = custom_data_bm.X;

for nn=1:length(p_inds)
    % Figure Name
    FigName = Y(p_inds(nn)).FigName;        

    
    % Name of data
    names{nn} = Y(p_inds(nn)).YName;
    
    % Assign Options
    bm_custom_opts.Name = Y(p_inds(nn)).YName;
    bm_custom_opts.FigLabel = FigLabel;
    bm_custom_opts.FigName = FigName;
    bm_custom_opts.FitFlags = FitFlags;
    bm_custom_opts.xstr = custom_data_bm.XStr;
    bm_custom_opts.Ind = nn;
    
    [hFs(nn),fouts{nn}] = customFit(X,Y(p_inds(nn)).Y,bm_custom_opts); 
%     ylim([0 10e4]);
    
    if ixon_doSave;saveFigure(hFs(nn),FigName,saveOpts);end

end

%% Save Data
    if ixon_doSave
        save([saveDir filesep 'custom_data_ixon_bm'],'custom_data_bm','fouts','names');
    end
    
    % if doSave && doUpload && exist(GDrive_root,'dir')
    %     gDir = [fileparts(getImageDir2(datevec(now),GDrive_root)) filesep FigLabel];
    %     gFile = [gDir filesep 'custom_data_ixon_bm'];        
    %     if ~exist(gDir,'dir')
    %        mkdir(gDir) 
    %     end
    %     save(gFile,'data');
    % end
