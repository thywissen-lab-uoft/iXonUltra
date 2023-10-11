% ixon_main_super.m
%
% Author : CF Fujiwara
%

% Display this filename
disp(repmat('-',1,60));disp(repmat('-',1,60));    
disp(['Calling ' mfilename '.m']);
disp(repmat('-',1,60));disp(repmat('-',1,60));    

% Add all subdirectories for this m file
curpath = fileparts(mfilename('fullpath'));
addpath(curpath);addpath(genpath(curpath));

a = fileparts(curpath);
addpath(a);addpath(genpath(a));

src = 'X:\Data';

%%

runs = [2023 10 07 33;
    2023 10 07 34;
    2023 10 07 35;
    2023 10 07 36;
    2023 10 07 37;
    2023 10 07 38;
    2023 10 07 39;
    2023 10 07 40;
    2023 10 07 41;
    2023 10 07 42;
    2023 10 07 43;
    2023 10 07 44;
    2023 10 07 45;
    2023 10 07 46;
    2023 10 07 47;
    2023 10 07 48;
    2023 10 07 49;
    2023 10 07 50;
    2023 10 07 51;
    2023 10 07 52;
    ];
%%
runs = [2023 10 09 03;
    2023 10 09 04;
2023 10 09 05;
2023 10 09 06;
2023 10 09 07;
2023 10 09 08;
2023 10 09 09;
2023 10 09 10;
2023 10 09 11;
2023 10 09 12;
2023 10 09 13;
2023 10 09 14;
2023 10 09 15;
    ];
%% Get the direcotry list
dir_list = ixon_findRunDirectory(runs);

for nn=1:length(dir_list)
   ixon_auto_dir = 0;
   imgdir = dir_list{nn};
   ixon_main;
    ixon_auto_dir = 1;
end



