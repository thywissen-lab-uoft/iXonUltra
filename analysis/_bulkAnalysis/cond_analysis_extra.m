runs = [
  2024 02 20 06;
  2024 02 20 07;
  2024 02 20 08;
  2024 02 20 09;
  2024 02 20 10;
  2024 02 20 11;
  2024 02 20 12;
  2024 02 20 13;
  2024 02 20 14;
  2024 02 20 15;
  2024 02 21 03;
  2024 02 21 04;


    ];

% runs = [
%   2024 02 15 05;
%   2024 02 15 06;
%   2024 02 15 07;
%   2024 02 15 08;
%   2024 02 15 09;
%   2024 02 15 10;
%   2024 02 15 11;
%     ];

% runs = [
%   2024 02 15 12;
%   2024 02 15 13;
%   2024 02 15 14;
%   2024 02 15 15;
%   2024 02 15 16;
%   2024 02 15 17;
%   2024 02 15 18;
%     ];

dir_list = ixon_findRunDirectory(runs);
files = {};
xdata = [];
ydata = [];
yerr = [];
for nn=1:length(dir_list)
   ixon_auto_dir = 0;
   imgdir = dir_list{nn};
   disp(['Loading data from ' imgdir]);
   files{nn}=load([imgdir filesep 'Figures' filesep 'conductivity_data.mat']);

   xdata(nn) = files{nn}.Params(1).conductivity_mod_ramp_time;
   ydata(nn) = abs(files{nn}.A);
   yerr(nn) = abs(files{nn}.Aerr);

end


f1 = figure(1005);
errorbar(xdata,ydata,yerr,'ko','MarkerFaceColor', 'k', 'DisplayName','4V');

xlabel('conductivity mod ramp time (ms)')
ylabel('amplitude response (sites)')
title('2V drive, 70 Hz, 190 G')



f2 = figure(1006);
errorbar(1./xdata,ydata,yerr,'ko','MarkerFaceColor', 'k', 'DisplayName','4V');

xlabel('conductivity mod ramp rate (kHz)')
ylabel('amplitude response (sites)')
title('2V drive, 70 Hz, 190 G')



