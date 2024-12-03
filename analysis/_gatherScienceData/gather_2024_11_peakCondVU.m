%THIS IS OBSOLETE: WE USE gather_2024_11_09_51Hz_PeakCondVersusU 

clear composite_data
composite_data = struct;
index=1;

%% Saving and Output

output_filename = '2024_11_peakCondVU_2.5Er_65mW';
doUpload = 1;
GDrive_root =['G:\My Drive\Lattice Shared\SharedData\Conductivity_Saturated_23-24'];

%% Define the Runs

%vary field 50 ms mod ramp 11/11 evap to 65 mW
composite_data(index).Name = '11/09 11/10 11/11';
composite_data(index).Description = '2.5Er, 50 ms mod ramp, 51 Hz 0.4 V & 0.8 V, vary field';
composite_data(index).Runs= [
       2024 11 09 14; %0.4V
       2024 11 09 15; %0.4V
       2024 11 09 16; %0.4V
       2024 11 09 17; %0.4V
       2024 11 09 18; %0.4V
       2024 11 09 19; %0.4V
       2024 11 09 20; %0.4V
       2024 11 09 21; %0.4V
       2024 11 09 22; %0.4V
       2024 11 09 23; %0.4V
       2024 11 09 24; %0.8V
       2024 11 09 25; %0.8V
       2024 11 09 26; %0.8V
       2024 11 09 27; %0.8V
       2024 11 10 01; %0.8V
       2024 11 10 02; %0.8V
       2024 11 10 03; %0.8V
       2024 11 10 04; %0.8V
       2024 11 10 05; %0.8V
       2024 11 10 06; %0.8V
       2024 11 10 07; %0.4V
       2024 11 10 08; %0.4V
       2024 11 10 09; %0.4V
       2024 11 10 10; %0.4V
       2024 11 10 11; %0.4V  
       2024 11 10 12; %0.4V
       2024 11 10 13; %0.4V
       2024 11 10 14; %0.4V
       2024 11 10 15; %0.4V
       2024 11 10 16; %0.4V
       2024 11 10 17; %0.8V
       2024 11 10 18; %0.8V
       2024 11 10 19; %0.8V
       2024 11 10 20; %0.8V
       2024 11 10 21; %0.8V
       2024 11 10 22; %0.8V
       2024 11 10 23; %0.8V
       2024 11 10 24; %0.8V
       2024 11 10 25; %0.8V
       2024 11 10 26; %0.8V
       2024 11 10 27; %0.4V
       2024 11 11 01; %0.4V 
       2024 11 11 02; %0.4V
%        2024 11 11 03; %0.4V - lost atoms halfway through
       2024 11 11 04; %0.2V
       2024 11 11 05; %0.2V
       2024 11 11 06; %0.2V
       2024 11 11 07; %0.2V 150 ms ramp
       2024 11 11 08; %0.2V 150 ms ramp
       2024 11 11 09; %0.2V 150 ms ramp
       2024 11 12 07; %0.2V
       2024 11 12 08; %0.2V
       2024 11 12 09; %0.2V
       2024 11 12 10; %0.2V 150 ms ramp
       2024 11 12 11; %0.2V 150 ms ramp
       
    ];
index=index+1;

%% Redo Analysis
do_redo_analysis = 0;    % Do you want to run analysis on it?

if do_redo_analysis
    opts=struct;
    opts.do_ixon_main           = 0;   % ixon_main
    opts.do_ixon_bin_analysis   = 0;   % ixon_bing
    opts.do_ixon_dig_analysis   = 1;   % ixon_dig
    ixon_super(composite_data,opts)
end

%% Gather Data
composite_data = gatherCompositeData(composite_data);

%% Upload the data

 try
     if ~exist(GDrive_root,'dir')
     mkdir(GDrive_root);
     end
 end
 
 if  doUpload && exist(GDrive_root,'dir')   
     fprintf('upload to google drive ...');
    gFile = fullfile(GDrive_root,[output_filename '.mat']);
    save(gFile,'composite_data'); 
    disp('done!')
 end