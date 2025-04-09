clear composite_data
composite_data = struct;

%% 200 G (requires different experimental trap parameters in analysis)
% runs = [  
%     2025 01 27 03;
%     2025 01 27 04;
%     2025 01 27 05;
%     2025 01 27 06;
%     2025 01 27 07;
%     2025 01 27 08;
%     2025 01 27 09;
%     2025 01 27 10;
%     2025 01 28 01;
%     2025 01 28 02;
%     2025 01 28 03;
%     2025 01 28 04;
%     2025 01 28 05;
%     2025 01 28 06;
%     2025 01 28 07;
%     ];
% 
% composite_data(end+1).Runs = runs;
% composite_data(end).Name = '2025_01_27 200 G specta 65 mW, 6 Er pulse';
%% 190 G (requires different experimental trap parameters in analysis)
runs = [  
    2025 01 28 08;
    2025 01 28 09;
    2025 01 28 10;
%     2025 01 28 11; %out of focus
    2025 01 28 12;
    2025 01 28 13;
    2025 01 28 14;
    2025 01 28 15;
    2025 01 28 16;
    2025 01 28 17;
    2025 01 28 18;
    2025 01 28 19;
    2025 01 28 20;
    2025 01 28 21;
    ];

composite_data(end+1).Runs = runs;
composite_data(end).Name = '2025_01_27 190 G specta 64.5 mW, 5.5 Er pulse';

%% 200.9 G
runs = [  
    2025 03 12 16;
    2025 03 12 17;
    2025 03 12 18;
    2025 03 12 19;
    2025 03 12 20;
    2025 03 12 21;
    2025 03 13 01;
    2025 03 13 02;
    2025 03 13 03;
    2025 03 13 04;
    2025 03 13 06;
    2025 03 13 07;
    2025 03 13 08;
    2025 03 13 09;
    2025 03 13 10;
    2025 03 13 11;
    ];

composite_data(end+1).Runs = runs;
composite_data(end).Name = '2025_03_12 200.9 G specta 54.5 mW, 1 Er pulse';

%% 200.6 G
runs = [
    2025 03 15 15
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
    2025 03 16 01
    2025 03 16 02
    2025 03 16 03
    2025 03 16 04
    2025 03 16 05
    2025 03 16 06
    2025 03 16 07
    2025 03 16 08
    2025 03 16 09
    2025 03 16 10
    2025 03 16 11
    2025 03 16 12
    2025 03 16 13
    ];

composite_data(end+1).Runs = runs;
composite_data(end).Name = '2025_03_16 201.1 G spectrum 54 mW, 3 Er pulse';

%% 198.5 G
runs = [
    2025 03 16 15
    2025 03 16 16
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
    2025 03 30 04 % 
    2025 03 30 05 %
    2025 03 30 06
    2025 03 30 07 %
    2025 03 30 08
    2025 03 30 09 %
    2025 03 30 10 %
    2025 03 30 11 %
    2025 03 30 12
    2025 03 30 13
    2025 03 30 14
    2025 03 30 15
    2025 03 30 16 %
    2025 03 30 17
    2025 03 30 18
    2025 03 31 01
    2025 03 31 02
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

%%
runs = [
    2025 04 08 20
    2025 04 08 21
    2025 04 08 22
    2025 04 08 23
    2025 04 08 24    
    2025 04 09 01
    2025 04 09 02
    2025 04 09 03 
    2025 04 09 04
    2025 04 09 05
    2025 04 09 06
    2025 04 09 07
    2025 04 09 08
    2025 04 09 09
    2025 04 09 10
    2025 04 09 11
    ];
composite_data(end+1).Runs = runs;
composite_data(end).Name = '2025_04_08 2.5 Er 200 G spectrum 53.7 mW, 2.5 Er pulse';

%% Gather All Data
composite_data(1)=[];
composite_data = gatherCompositeData(composite_data);
