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
composite_data(index).Type = 'peak';
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
        2024 12 03 15;

    ];
index=index+1;
%% 2024/12/02-2024/12/03

%Single frequency 54 Hz, 0.8 V, Variable Field
% 2.5 Er, 50 ms mod ramp, 68 mW
composite_data(index).Name = '12/02 0.8V 68 mW 54 Hz';
composite_data(index).Description = '12/02 0.8V 68 mW 54 Hz';
composite_data(index).Type = 'peak';
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
        2024 12 03 14;
        2024 12 03 16;
    ];
index=index+1;

%% 2024/12/03-2024/12/04
% 
% %Single frequency 54 Hz, vary field evap depth and pulse depth
% 
% % 2.5 Er, 50 ms mod ramp, 190 G, 64-70 mW, 0 Er
% composite_data(index).Name = '2024_12_03 54 Hz, 190 G, 64-70 mW, 0 Er, 0.4 V';
% composite_data(index).Description = '54 Hz, 2.5 Er, 50 ms mod ramp, 190 G, 64-70 mW, 0 Er pulse, 0.4 V';
% composite_data(index).Type = 'peak';
% composite_data(index).Runs = [     
%     2024 12 03 19;
%     2024 12 03 21;
%     2024 12 03 23;
%     2024 12 04 02;
%     2024 12 04 12;
%     ];
% 
% index=index+1;
% 
% % 2.5 Er, 50 ms mod ramp, 190 G, 64-70 mW, 4 Er pulse
% composite_data(index).Name = '2024_12_03 54 Hz, 190 G, 64-70 mW, 4 Er, 0.4 V';
% composite_data(index).Description = '54 Hz, 2.5 Er, 50 ms mod ramp, 190 G, 64-70 mW, 4 Er pulse, 0.4 V';
% composite_data(index).Type = 'peak';
% composite_data(index).Runs = [     
%     2024 12 03 18;
%     2024 12 03 20;
%     2024 12 03 22;
%     2024 12 04 01;
%     2024 12 04 11;
%     2024 12 04 13;
%     ];
% 
% index=index+1;
% 
% % 2.5 Er, 50 ms mod ramp, 200G, 64-70 mW, 0 Er pulse
% composite_data(index).Name = '2024_12_03 54 Hz, 200 G, 64-70 mW, 0 Er, 0.8 V';
% composite_data(index).Description = '54 Hz, 2.5 Er, 50 ms mod ramp, 200G, 64-70 mW, 0 Er pulse, 0.8 V';
% composite_data(index).Type = 'peak';
% composite_data(index).Runs = [     
%     2024 12 04 04;
%     2024 12 04 06;
%     2024 12 04 08;
%     2024 12 04 10;
%     ];
% 
% index=index+1;
% 
% % 2.5 Er, 50 ms mod ramp, 200G, 64-70 mW, 4 Er pulse
% composite_data(index).Name = '2024_12_03 54 Hz, 200 G, 64-70 mW, 4 Er, 0.8 V';
% composite_data(index).Description = '54 Hz, 2.5 Er, 50 ms mod ramp, 200G, 64-70 mW, 4 Er pulse, 0.8 V';
% composite_data(index).Type = 'peak';
% composite_data(index).Runs = [     
%     2024 12 04 03;
%     2024 12 04 05;
%     2024 12 04 07;
%     2024 12 04 09;
%     ];
% 
% index=index+1;

%% Redo Analysis
do_redo_analysis =0;    % Do you want to run analysis on it?

if do_redo_analysis
    opts=struct;
    opts.do_ixon_main           = 0;   % ixon_main
    opts.do_ixon_bin_analysis   = 0;   % ixon_bing
    opts.do_ixon_dig_analysis   = 1;   % ixon_dig
    ixon_super(composite_data,opts)
end

%% Gather Data
composite_data = gatherCompositeData(composite_data);

