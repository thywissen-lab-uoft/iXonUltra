clear composite_data
composite_data = struct;
index=1;
%% Introduction

%% 11/26-11/27
% composite_data= struct;
% composite_data.Name = '2024_11_26 201.1 G 1.5 V Drive, Full Spectra, 70 mW Evap';
% composite_data.Runs = [     
%     2024 11 26 12;
%     2024 11 26 13;
%     2024 11 26 14;
%     2024 11 26 16;
%     2024 11 26 16;
%     2024 11 26 17;
%     2024 11 26 18;
%     2024 11 26 19;
%     2024 11 27 01;
%     2024 11 27 02;
%     2024 11 27 03;
%     2024 11 27 04;
%     ];
% index=index+1;

%% 2024/11/28-2024/11/29
%201.1 G high field 50 ms mod ramp 11/22-11/23 vary evap depth
composite_data(index).Name = '11/28-11/29 201.1 G 4 V';
composite_data(index).Description = '201.1 G high field, 2.5Er, 50 ms mod ramp, 4V, Full Spectrum';
composite_data(index).Runs= [ 
        2024 11 28 04;
        2024 11 28 05;
        2024 11 28 06;
        2024 11 28 07;
        2024 11 28 08;
        2024 11 28 09;
        2024 11 28 10;
        2024 11 28 11;        
        2024 11 28 12;        
        2024 11 29 01;
        2024 11 29 02;
        2024 11 29 03;
        2024 11 29 04;
        2024 11 29 05;
        2024 11 29 06;
        2024 11 29 07;
        2024 11 29 08;
        2024 11 29 09;
        2024 11 29 10;
                
    ];
index=index+1;

%%

composite_data(index).Name = '12_01-02 201.1 G 4 V,full-spec,66.5 mW,5 Er pulse';
composite_data(index).Runs =[     
    2024 11 29 19;
    2024 11 29 20;
    2024 11 29 21;
    2024 11 29 22;
    2024 11 29 23;
    2024 11 29 24;
    2024 11 29 25;
    2024 11 29 26;
    2024 11 30 01;
    2024 11 30 02;
    2024 11 30 03;
    2024 11 30 04;
    2024 11 30 05;
    2024 11 30 06;
    2024 11 30 07;
    2024 11 30 08;
    2024 11 30 09;
    2024 11 30 10;
    2024 11 30 11;
    2024 11 30 12;
    2024 11 30 13;
    2024 11 30 14;
    2024 11 30 15;
    2024 11 30 16;
    2024 11 30 17;
    2024 11 30 18;

    ];
index=index+1;


%% 2024/12/01-2024/12/02
composite_data(index).Name = '2024_12_01-02 201.1 G 4 V Drive, Full Spectra, 67.5 mW Evap, 5 ER Pulse';
composite_data(index).Runs = [     
    2024 12 01 02;
    2024 12 01 03;
    2024 12 01 04;
    2024 12 01 05;
    2024 12 01 06;
    2024 12 01 07;
    2024 12 01 08;
    2024 12 01 09;
    2024 12 01 10;
    2024 12 01 11;
    2024 12 01 12;
    2024 12 01 13;
    2024 12 01 14;
    2024 12 01 15;
    2024 12 01 16;
    2024 12 02 01;
    2024 12 02 02;
    2024 12 02 03;
    2024 12 02 04;
    2024 12 02 05;
    2024 12 02 06;
    2024 12 02 07;
    2024 12 02 08;
    2024 12 02 09;
    2024 12 02 10;
    2024 12 02 11;
    2024 12 02 12;
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

