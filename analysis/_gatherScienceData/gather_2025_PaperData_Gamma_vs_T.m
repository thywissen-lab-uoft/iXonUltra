clear composite_data
composite_data = struct;
warning off
%% T = 1t (201.1G) 
runs = [  
    2024 11 08 19;
    2024 11 08 20;
    2024 11 08 21;
    2024 11 08 22;
    2024 11 08 23;
    2024 11 08 24;
    2024 11 08 25;
    2024 11 08 26;
    2024 11 09 01;
    2024 11 09 02; % small signal
    2024 11 09 03;
    2024 11 09 04;
    2024 11 09 05;
    2024 11 09 06;
    2024 11 09 07;
    2024 11 09 08;
    2024 11 09 09;
    2024 11 09 10;
    2024 11 09 11;
    ];

composite_data(end+1).Runs = runs;
composite_data(end).Name = '2024_11_08 201.1 G spectra T = 1t';
%% T = 1t (201.1G) 
runs = [  
%     2024 12 05 01; % plane selection issues? inconsistent density
    2024 12 05 02; % increasing atom number and density, but overall OK
    2024 12 05 03; % OK, fluctuating atom number
    2024 12 05 04; % OK, a few bad shots
    2024 12 05 05; % OK
    2024 12 06 01; % OK
    2024 12 06 02; % inconsistent/low density
    2024 12 06 03; % OK
    2024 12 06 04; % OK
    2024 12 06 05; % OK
    2024 12 06 06; % a few shots of low density/atom number
    2024 12 06 07; % inconsistent/low density
    2024 12 06 08; % OK
    2024 12 06 09; % decreasing density
    2024 12 06 10; % OK
%     2024 12 06 11; % inconsistent/low density & low amplitude response
    2024 12 06 12; 
    ];

composite_data(end+1).Runs = runs;
composite_data(end).Name = '2024_12_05 201.1 G spectra T = 1t';
%% T = 2.9t (201.1G) 
runs = [  
    2024 12 09 02; % OK
    2024 12 09 03; % OK
    2024 12 09 04; % fluctuating density, low amplitude response/bad fit
%     2024 12 09 05; % decreasing density, slightly out of focus, bad fit
    2024 12 09 06; % not great signal, but OK
%     2024 12 10 01; % plane selection issues? inconsistent density, bad fit
    2024 12 10 02; % not great fit
    2024 12 10 03; % OK
    2024 12 10 04; % low density, a few bad shots
    2024 12 10 05; % OK
    2024 12 10 06; % OK
    2024 12 10 07; % OK
    2024 12 10 08; % OK
    2024 12 10 09; % OK
%     2024 12 10 10; % bad fit, inconsistent density
    2024 12 10 11; % OK
    2024 12 10 12; % OK
    2024 12 10 13; % OK
    2024 12 10 14; % OK
    ];

composite_data(end+1).Runs = runs;
composite_data(end).Name = '2024_12_09 201.1 G spectra T = 2.9t';

%% T = 1.7t (201.1G) 
runs = [  
    2025 01 23 07;
    2025 01 23 08;
    2025 01 23 09;
    2025 01 23 10;
    2025 01 23 11;
    2025 01 23 12;
    2025 01 24 01;
    2025 01 24 02;
%     2025 01 24 03; % no signal
    2025 01 24 04;
    2025 01 24 05;
%     2025 01 24 06; % no signal, low density
    2025 01 24 07; % OK, but slightly out of focus
    2025 01 24 08; % not the best fit, slightly out of focus
    2025 01 24 09;
    2025 01 24 10;
    2025 01 24 11;
    2025 01 24 12;
    ];

composite_data(end+1).Runs = runs;
composite_data(end).Name = '2025_01_23 201.1 G spectra T = 1.7t';

%% T = 0.93t (201.1G) 
runs = [  
    2025 01 25 25;
    2025 01 25 26;
%     2025 01 25 27; % no signal
    2025 01 26 01; 
    2025 01 26 02; 
    2025 01 26 03;
    2025 01 26 04;
    2025 01 26 05; % not the best signal, slightly higher density
    2025 01 26 06; 
    2025 01 26 07;
    2025 01 26 08; % this is good, not sure why it was removed previously
    2025 01 26 09; 
    2025 01 26 10;
    2025 01 26 11;
    2025 01 26 12;
    2025 01 26 13;
    2025 01 26 14;
    2025 01 26 15;
    2025 01 26 16;
    2025 01 26 17;
    2025 01 26 18;
    2025 01 26 19;
    2025 01 26 20; % not great signal, fluctuating density
%     2025 01 26 21; % out of focus, binning & digitization is wrong
%     2025 01 26 22; % out of focus, binning & digitization is wrong
%     2025 01 26 23; % out of focus, binning & digitization is wrong
    2025 01 26 24; % back in focus
    2025 01 26 25; 
    2025 01 26 26;
    ];

composite_data(end+1).Runs = runs;
composite_data(end).Name = '2025_01_25 201.1 G spectra T = 0.93t';

%% T = 2.1t (201.1G) 
runs = [  
    2025 03 15 31;
    2025 03 15 32;
    2025 03 16 01; % decreasing density, no the best signal
    2025 03 16 02; 
    2025 03 16 03;
    2025 03 16 04;
    2025 03 16 05; % low density
    2025 03 16 06; 
    2025 03 16 07;
    2025 03 16 08; % not the best signal, heating, decreasing density
%     2025 03 16 09; % many bad shots
    2025 03 16 10;
    2025 03 16 11;
    2025 03 16 12;
    2025 03 16 13; % low density
    ];

composite_data(end+1).Runs = runs;
composite_data(end).Name = '2025_03_15 201.1 G spectra T = 0.93t';

%% T = 1.4t (201.1G) 
runs = [  
    2025 03 18 04;
    2025 03 18 05;
    2025 03 18 06;
    2025 03 18 07;
    2025 03 18 08;
    2025 03 18 09;
    2025 03 18 10;
    2025 03 18 11;
    2025 03 18 12;
    2025 03 18 13;
    2025 03 18 14; % not the best signal, but probably OK 
    2025 03 18 15; % a few bad shots
    2025 03 18 16;
    2025 03 18 17;
    2025 03 18 18;
    ];

composite_data(end+1).Runs = runs;
composite_data(end).Name = '2025_03_18 201.1 G spectra T = 1.4t';

%% T = 2.2t (201.1G) 
runs = [  
    2025 03 20 06; % density decreases, but good signal
    2025 03 20 07; % one bad shot
    2025 03 20 08; % fluctuating density, not best signal
    2025 03 20 09; % OK, but starting to go out of focus
    2025 03 21 01;
    2025 03 21 02;
    2025 03 21 03;
    2025 03 21 04;
    2025 03 21 05;
    2025 03 21 06;
    2025 03 21 07;
    2025 03 21 08;
    2025 03 21 09;
    2025 03 21 10;
    2025 03 21 11;
    2025 03 21 12;
    2025 03 21 13;
    2025 03 21 14;
    2025 03 21 15;
    2025 03 21 17;
    2025 03 21 22;
    2025 03 21 23;
    ];

composite_data(end+1).Runs = runs;
composite_data(end).Name = '2025_03_20 201.1 G spectra T = 2.2t';
%% Gather All Data
composite_data(1)=[];
composite_data = gatherCompositeData(composite_data);

%% Upload

% doUpload = true;
doUpload = false;

GDrive_root =['G:\My Drive\Lattice Shared\SharedData\Conductivity_Saturated_23-24'];
output_folder_name ='2025_paper_data_vs_T';
saveDir = fullfile(GDrive_root,output_folder_name);

if doUpload
    try
        if ~exist(GDrive_root,'dir');mkdir(GDrive_root);end
        if ~exist(saveDir,'dir');mkdir(saveDir);end
         gFile = fullfile(saveDir,'composite_data.mat');

        disp(gFile);
        fprintf('upload to google drive ...');        
        save(gFile,'composite_data'); 
        disp('done!')   
    catch ME
        error('OH NO I COUD NOT UPL<OD');
    end
end 