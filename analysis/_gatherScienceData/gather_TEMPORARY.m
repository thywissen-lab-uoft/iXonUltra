clear composite_data
composite_data = struct;
index=1;
%% Introduction


% %% 2024/11/28-2024/11/29
% %201.1 G high field 50 ms mod ramp 11/22-11/23 vary evap depth
% composite_data(index).Name = '11/28-11/29 201.1 G 4 V';
% composite_data(index).Description = '201.1 G high field, 2.5Er, 50 ms mod ramp, 4V, Full Spectrum';
% composite_data(index).Runs= [ 
%         2024 11 28 04;
%         2024 11 28 05;
%         2024 11 28 06;
%         2024 11 28 07;
%         2024 11 28 08;
%         2024 11 28 09;
%         2024 11 28 10;
%         2024 11 28 11;        
%         2024 11 28 12;
%         
%         2024 11 29 01;
%         2024 11 29 02;
%         2024 11 29 03;
%         2024 11 29 04;
%         2024 11 29 05;
%         2024 11 29 06;
%         2024 11 29 07;
%         2024 11 29 08;
%         2024 11 29 09;
%         2024 11 29 10;
%                 
%     ];
% index=index+1;

%% 2024/12/02-2024/12/03

%Single frequency 54 Hz, 0.4 V, Variable Field
% 2.5 Er, 50 ms mod ramp, 68 mW
composite_data(index).Name = '12/02 0.4V 68 mW 54 Hz';
composite_data(index).Description = '12/02 0.4V 68 mW 54 Hz';
composite_data(index).Runs= [ 
        2024 12 02 16;
        2024 12 02 18
        2024 12 02 20;
        2024 12 03 01;
        2024 12 03 03;
        2024 12 03 05;
        2024 12 03 07;
        2024 12 03 09;
        2024 12 03 11;
        2024 12 03 13;
    ];
index=index+1;
%% 2024/12/02-2024/12/03

%Single frequency 54 Hz, 0.8 V, Variable Field
% 2.5 Er, 50 ms mod ramp, 68 mW
composite_data(index).Name = '12/02 0.8V 68 mW 54 Hz';
composite_data(index).Description = '12/02 0.8V 68 mW 54 Hz';
composite_data(index).Runs= [ 
        2024 12 02 17;
        2024 12 02 19
        2024 12 02 21;
        2024 12 03 02;
        2024 12 03 04;
        2024 12 03 06;
        2024 12 03 08;
        2024 12 03 10;
        2024 12 03 12;
    ];
index=index+1;


%% Redo Analysis
do_redo_analysis = 0;    % Do you want to run analysis on it?

if do_redo_analysis
    opts=struct;
    opts.do_ixon_main           = 1;   % ixon_main
    opts.do_ixon_bin_analysis   = 1;   % ixon_bing
    opts.do_ixon_dig_analysis   = 1;   % ixon_dig
    ixon_super(composite_data,opts)
end

%% Gather Data
composite_data = gatherCompositeData(composite_data);

