clear composite_data
composite_data = struct;

%% 190 G (requires different experimental trap parameters in analysis)
runs = [  
    2025 01 28 08; % not great signal
    % 2025 01 28 09; % no signal
    2025 01 28 10;
%     2025 01 28 11; %out of focus
    2025 01 28 12;
    2025 01 28 13; % a bit larger
    2025 01 28 14;
    % 2025 01 28 15; % pretty much no signal, slightly out of focus
    2025 01 28 16;
    2025 01 28 17; % density a bit lower
    2025 01 28 18; % not much of a signal, slightly out of focus
%     2025 01 28 19; % no signal
    2025 01 28 20; % OK, a little warmer
    2025 01 28 21;
    ];

composite_data(end+1).Runs = runs;
composite_data(end).Name = '2025_01_28 190 G spectrum 64.5 mW, 5.5 Er pulse';

%% 200.9 G
runs = [  
    2025 03 12 16;
    2025 03 12 17;
    2025 03 12 18;
    2025 03 12 19;
    2025 03 12 20;
    2025 03 12 21;
    2025 03 13 01; % fluctuating density, inconsistent fluor/atom
    2025 03 13 02;
    2025 03 13 03; % decreasing density
    2025 03 13 04; % fluctuating density
    2025 03 13 05; % perfectly fine, not sure why this wasn't analyzed
    2025 03 13 06; 
    2025 03 13 07;
    2025 03 13 08;
    2025 03 13 09; % drop in density
    2025 03 13 10; 
    2025 03 13 11;
    ];

composite_data(end+1).Runs = runs;
composite_data(end).Name = '2025_03_12 200.9 G specta 54.5 mW, 1 Er pulse';

%% 200.6 G
runs = [
    2025 03 15 15 % decreasing density
    2025 03 15 16 
    2025 03 15 17
    2025 03 15 18
    2025 03 15 19 
    2025 03 15 20
    2025 03 15 21
    2025 03 15 22
    2025 03 15 23
    2025 03 15 24
    2025 03 15 25
    2025 03 15 26
    2025 03 15 27
    2025 03 15 28
    2025 03 15 29
    ];

composite_data(end+1).Runs = runs;
composite_data(end).Name = '2025_03_15 200.6 G spectrum 54 mW, 2.5 Er pulse';

%% 201.1 G
runs = [
    2025 03 15 31 
    2025 03 15 32
    2025 03 16 01 % 70 Hz density a little low
    2025 03 16 02
    2025 03 16 03
    2025 03 16 04
    2025 03 16 05 % 100 Hz lower density
    2025 03 16 06 
    2025 03 16 07
    2025 03 16 08 % not the best signal
%     2025 03 16 09 % few bad shots, no signal
    2025 03 16 10
    2025 03 16 11
    2025 03 16 12
    2025 03 16 13 % fluctuating/low density
    ];

composite_data(end+1).Runs = runs;
composite_data(end).Name = '2025_03_16 201.1 G spectrum 54 mW, 3 Er pulse';

%% 198.5 G
runs = [
    2025 03 16 15
%     2025 03 16 16 % large drift/no signal
    2025 03 16 17 
    2025 03 16 18
    2025 03 16 19
    2025 03 16 20
    2025 03 16 21
    2025 03 16 22
    2025 03 16 23
    2025 03 16 24
    2025 03 16 25
    2025 03 16 26
    2025 03 16 27
    2025 03 16 28
    2025 03 16 29
    ];

composite_data(end+1).Runs = runs;
composite_data(end).Name = '2025_03_16 198.5 G spectrum 54.3 mW, 3 Er pulse';

%% 195 G
runs = [
%     2025 03 30 04 % few bad shots, no signal
    2025 03 30 05 % few bad shots, OK signal
    2025 03 30 06 % not great signal, but still usable
    2025 03 30 07 % few bad shots, OK signal 
    2025 03 30 08
    2025 03 30 09 % low amplitude, but still usable
    2025 03 30 10 % decreasing density
%     2025 03 30 11 % no signal
    2025 03 30 12
    2025 03 30 13
    2025 03 30 14 % decreasing density
    2025 03 30 15 % fluctuating density
    2025 03 30 16 % decreasing density
    2025 03 30 17 
    2025 03 30 18 % lower density
    2025 03 31 01 % decreasing density
    2025 03 31 02 % decreasing density
    2025 03 31 03
    2025 03 31 05 %repeat
    2025 03 31 06 %repeat
    2025 03 31 07 %repeat
    2025 03 31 08 %repeat
    2025 03 31 09 %repeat
    2025 03 31 10 %repeat
    2025 03 31 11 %repeat
    2025 03 31 12 %repeat
    ];

composite_data(end+1).Runs = runs;
composite_data(end).Name = '2025_03_30 195 G spectrum 54 mW, 3 Er pulse';

%% 200 G 70 Hz trap frequency
runs = [
    2025 04 10 17 % okay signal
    2025 04 10 18 % large amplitude
    2025 04 10 19 % okay signal
    2025 04 10 20 % good signal
    2025 04 10 21 % good signal, large amplitude
    2025 04 10 22 % good signal, large amplitude
    2025 04 11 01 % good signal, large amplitude
    2025 04 11 02 % good signal
    2025 04 11 03 % good signal, large amplitude
    2025 04 11 04 % good signal
    2025 04 11 05 % good signal, large amplitude
    2025 04 11 06 % good signal, large amplitude
    2025 04 11 07 % good signal, small amplitude
    2025 04 11 08 % good signal
    2025 04 11 09 % good signal
%     2025 04 11 10 % not great signal, large drift
    ];
composite_data(end+1).Runs = runs;
composite_data(end).Name = '2025_04_10 2.5 Er 200 G spectrum 53.5 mW, 2.8 Er pulse 70 Hz';

%% Gather All Data
composite_data(1)=[];
composite_data = gatherCompositeData(composite_data);

%% Upload

doUpload = true;
% doUpload = false;

GDrive_root =['G:\My Drive\Lattice Shared\SharedData\Conductivity_Saturated_23-24'];
output_folder_name ='2025_paper_data_vs_U';
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