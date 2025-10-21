%% Export digdata for collabotors


%% Select image directory
    
if ~exist('dig_auto_file');dig_auto_file = 1;end

if dig_auto_file
    dialog_title='Select digdata';       
    ixon_getYearDir = 'X:\Data\2025\';
    [filename,dig_imgdir,b]=uigetfile(fullfile(ixon_getYearDir,'*.mat'),dialog_title);
    filename = fullfile(dig_imgdir,filename);    
    if  b == 0
        disp('Canceling.');    
        return; 
    end
end
% Load the Data
clear digdata
tic
fprintf('Loading digdata ...');
digdata = load(filename);
disp([' done (' num2str(toc,'%.2f') 's)']);

disp(mean(digdata.Ys_site./digdata.Xs_site));

%% Shift all data by the COM

fprintf('Assigning to output...\n');

iCOM = round(digdata.Xc_site,0)-100; % shift to site 100

for i=1:length(digdata.Xc_site)
    
    digdata.ZdigShift(:,:,i) = circshift(digdata.Zdig(:,:,i),-iCOM(i),2);
end

% increase Zdig by fidelity
fidelity = 0.88;
digdata.ZdigShift = digdata.ZdigShift/fidelity;

% Make array of average data

out = struct;

% Occupation
out.Nmean = mean(digdata.ZdigShift,3);
out.Nstd  = std(digdata.ZdigShift,0,3);

% trap frequency (taken from lineshape)
out.trap_freq = 67;
out.trap_freq_unc = 3;
out.trap_freq_units = 'Hz';

% Interaction strength
out.U = 3229; % for 201.1 G + 0.1238 G
% out.U = 1683; % for 200 G + 0.1238 G
% out.U = 757; % for 195 G + 0.1238 G
% out.U = 1130; % for 198.5 + 0.1238 G
% out.U = 2231; % for 200.6 + 0.1238 G
% out.U = 2723; % for 200.9 + 0.1238 G
% out.U = 601; % 190 G + 0.1238 G
out.U_units = 'Hz';

% tunneling rate
out.t = 568; % for VL = 2.47 Er
out.t_unc = 11; % full range
out.t_units = 'Hz';

% constants
out.m = 6.6422e-26;
out.m_units = 'kg';
out.a_L = 527;
out.a_L_units = 'nm';

%% Save output

saveDir = 'Y:\Data for colab\';
filename = 'Data_201_1G_03_20_80Hz.mat';

try if ~exist(saveDir,'dir');mkdir(saveDir);end;end
    filename = fullfile(saveDir,filename);
    disp(['Saving ' filename ' ...']);
    save(filename, '-struct','out');

