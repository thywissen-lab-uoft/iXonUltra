function [ixondata,qpd_out] = AnalyzeIxonQPD(ixondata,saveDir,opts)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matches QPD files to a list of iXon data of the form ixondata.digdata
% Analyzes each run to get average powers, fit parameters for data list

% opts fields are:
% doPlot (plots QPD analyis)
% doSave (saves plots to iXon folder)
% isMac (sets source dir.)
% isRemote 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Display

disp('Analyzing QPD traces...')

%% Initialize options

if nargin~=3
    opts.doPlot = 1;
    opts.isMac = 0;
    opts.doSave = 1;
    opts.isRemote = 0;
end

%% Close plots
%  CF removing this because you dont know the figure number a priori
% if ishandle(101)
%     close(101)
% end
% if ishandle(201)
%     close(201)
% end
% if ishandle(301)
%     close(301)
% end
% if ishandle(401)
%     close(401)
% end

%% Data Root
% Find the data source file

if opts.isMac       
    pdsrc =  '/Volumes/main/LabJackLogs/ODTQPD'; % Mac
else
    pdsrc = 'X:\LabJackLogs\ODTQPD'; %Windows
end


%% Load the QPD files
qpdfiles=strings(0);

for jj=1:length([ixondata.Params])

    % Date string for run
    fdate = datestr(ixondata(jj).Params.ExecutionDateStr,'YYYY-mm-dd_HH-MM-SS');
    
    % Directory of QPD file
    pddir = fullfile(pdsrc,fdate(1:4),[fdate(1:4) '.' fdate(6:7)],[fdate(6:7) '.' fdate(9:10)]);
    
    % QPD filename
    fname = ['ODTQPD_' fdate '.mat'];
    qpdfiles(end+1) = fullfile(pddir,fname);

end

%% Analyze the QPD traces

[ixondata,qpd_out] = loadQPDruns(ixondata,qpdfiles,opts);


%% Plotting

% Plot the individual fits during modulation and summary plot of modulation traces
if opts.doPlot

    for nn=1:length(qpdfiles)

   
        qpd_data(nn) = ixondata(nn).qpd_data;

        %%%%%%%%%%%%% Plot individual modulation traces %%%%%%%%%%%%%%%

        tt=linspace(min(qpd_data(nn).modfit_t),max(qpd_data(nn).modfit_t),1e3);
    
        f=figure(9001);
        f.Position = [100 100 540 400]; 
        clf;
    
        % Plot X1
        subplot(211)
        plot(qpd_data(nn).modfit_t,qpd_data(nn).modfit_normX1data)
        hold on;
        plot(tt,feval(qpd_data(nn).modfit_X1,tt),'-')
        hold off;
    
        xlim([0 50])
        xlabel('Time (ms)')
        ylabel('Normalized X1')
    
        % Plot X2
        subplot(212)
        plot(qpd_data(nn).modfit_t,qpd_data(nn).modfit_normX2data)
        hold on;
        plot(tt,feval(qpd_data(nn).modfit_X2,tt))
        hold off;
    
        xlim([0 50])
        xlabel('Time (ms)')
        ylabel('Normalized X2')
    

        %%%%%%%%%%%% Plot summary of modulation traces %%%%%%%%%%%%%%

        % Select the data during the modulation, from ramp up to end
%         ramp_start = find(round(qpd_data(nn).t,1)==-150);
        ramp_start = find(round(qpd_data(nn).t,1)==-ixondata(nn).Params.conductivity_mod_ramp_time);

        % Find end of mod from lattice ramp up BUGGY RL FIX
        % (approx. as max slope of lattice ramp)
%         mod_end   = find(diff(qpd_data(nn).XLATT)==max(diff(qpd_data(nn).XLATT)));
        mod_end   = max(find(diff(qpd_data(nn).TRIG)==-1));

        % Add 50 ms buffer time to each end 
        buffer_t = 250;
    
        fmsum=figure(201);
        fmsum.Position = [600 100 800 600];
        plotColor = jet(length(qpdfiles));

        ax1=subplot(4,3,1);
        plot(qpd_data(nn).t,qpd_data(nn).X1./qpd_data(nn).SUM1,'.',  'color',  plotColor(nn,:))
        hold on;
        xlabel('Time (ms)')
        ylabel('ODT 1 X (V/V)')
    
        ax2=subplot(4,3,4);
        plot(qpd_data(nn).t,qpd_data(nn).Y1./qpd_data(nn).SUM1,'.',  'color',  plotColor(nn,:))
        hold on;
        xlabel('Time (ms)')
        ylabel('ODT 1 Y (V/V)')
    
        ax3=subplot(4,3,7);
        plot(qpd_data(nn).t,qpd_data(nn).X2./qpd_data(nn).SUM2,'.',  'color',  plotColor(nn,:))
        hold on;
        xlabel('Time (ms)')
        ylabel('ODT 2 X (V/V)')
    
        ax4=subplot(4,3,10);
        plot(qpd_data(nn).t,qpd_data(nn).Y2./qpd_data(nn).SUM2,'.',  'color',  plotColor(nn,:))
        hold on;
        xlabel('Time (ms)')
        ylabel('ODT 2 Y (V/V)')
    
        linkaxes([ax1 ax2 ax3 ax4],'x')
        xlim([min(qpd_data(nn).t(ramp_start-buffer_t:mod_end+buffer_t)) max(qpd_data(nn).t(ramp_start-buffer_t:mod_end+buffer_t))])
        ylim('auto')    
        
        % Plot modulation in um

        fmsum_um=figure(401);
        fmsum_um.Position = [600 100 800 600];
        plotColor = jet(length(qpdfiles));
    
        ax1=subplot(4,3,1);
        plot(qpd_data(nn).t,qpd_data(nn).X1./qpd_data(nn).SUM1*94.0311,'.',  'color',  plotColor(nn,:))
        hold on;
        xlabel('Time (ms)')
        ylabel('ODT 1 X (um)')
    
        ax2=subplot(4,3,4);
        plot(qpd_data(nn).t,qpd_data(nn).Y1./qpd_data(nn).SUM1,'.',  'color',  plotColor(nn,:))
        hold on;
        xlabel('Time (ms)')
        ylabel('ODT 1 Y (V/V)')
    
        ax3=subplot(4,3,7);
        plot(qpd_data(nn).t,qpd_data(nn).X2./qpd_data(nn).SUM2*994.9276,'.',  'color',  plotColor(nn,:))
        hold on;
        xlabel('Time (ms)')
        ylabel('ODT 2 X (um)')
    
        ax4=subplot(4,3,10);
        plot(qpd_data(nn).t,qpd_data(nn).Y2./qpd_data(nn).SUM2,'.',  'color',  plotColor(nn,:))
        hold on;
        xlabel('Time (ms)')
        ylabel('ODT 2 Y (V/V)')
    
        linkaxes([ax1 ax2 ax3 ax4],'x')
        xlim([min(qpd_data(nn).t(ramp_start-buffer_t:mod_end+buffer_t)) max(qpd_data(nn).t(ramp_start-buffer_t:mod_end+buffer_t))])
        ylim('auto')    
   
    end

    %%%%%%%%%%%%%% Summary plot for the modulation %%%%%%%%%%%%%%%

    % (V/V)

    figure(201)
    ax5=subplot(4,3,[2 5]);
    histogram([qpd_data.modfit_X1_A],15)
    xlabel('Norm. X1 mod. amp. (V/V)')
    
    ax6=subplot(4,3,[3 6]);
    histogram([qpd_data.modfit_X2_A],15)
    xlabel('Norm. X2 mod. amp. (V/V)')
    
    ax7=subplot(4,3,[8 11]);
    histogram([qpd_data.modfit_X1_B]/(2*pi),15)
    hold on;
    histogram([qpd_data.modfit_X2_B]/(2*pi),15)
    xlabel('Phase/ (2$\pi$)','Interpreter','Latex','FontName','Helvetica')
    
    ax8=subplot(4,3,[9 12]);
    histogram(1./[qpd_data.modfit_X1_T],15)
    hold on;
    histogram(1./[qpd_data.modfit_X2_T],15)
    xlabel('Frequency (kHz)')

    % um
    figure(401)

    % um Feb 13/24 conversion 94.0311 um/(V/V)
    ax5=subplot(4,3,[2 5]);
    histogram([qpd_data.modfit_X1_A]*94.0311,15)
    xlabel('Norm. X1 mod. amp. (um)')
    
    % um Feb 13/24 conversion 994.9276 um/(V/V)
    ax6=subplot(4,3,[3 6]);
    histogram([qpd_data.modfit_X2_A]*994.9276,15)
    xlabel('Norm. X2 mod. amp. (um)')
    
    ax7=subplot(4,3,[8 11]);
    histogram([qpd_data.modfit_X1_B]/(2*pi),15)
    hold on;
    histogram([qpd_data.modfit_X2_B]/(2*pi),15)
    xlabel('Phase/ (2$\pi$)','Interpreter','Latex','FontName','Helvetica')
    
    ax8=subplot(4,3,[9 12]);
    histogram(1./[qpd_data.modfit_X1_T],15)
    hold on;
    histogram(1./[qpd_data.modfit_X2_T],15)
    xlabel('Frequency (kHz)')
    
    %%%%%%%%%%%%% Summary plot for the power %%%%%%%%%%%%%
    
    fpsum = figure(301);
    
    axODT1A=subplot(3,4,1);
    histogram([qpd_data.ODT1_ave],15)
    xlabel('Sum 1 Ave. (V)')
    
    axODT1SD=subplot(3,4,2);
    histogram([qpd_data.ODT1_std],15)
    xlabel('Sum 1 Std. Dev. (V)')
    
    axODT2A=subplot(3,4,3);
    histogram([qpd_data.ODT2_ave],15)
    xlabel('Sum 2 Ave. (V)')
    
    axODT2SD=subplot(3,4,4);
    histogram([qpd_data.ODT2_std],15)
    xlabel('Sum 2 Std. Dev. (V)')
    
    axLXA=subplot(3,4,5);
    histogram([qpd_data.XLatt_ave],15)
    xlabel('X Latt. Ave. (V)')
    
    axLXSD=subplot(3,4,6);
    histogram([qpd_data.XLatt_std],15)
    xlabel('X Latt. Std. Dev. (V)')
    
    axLYA=subplot(3,4,7);
    histogram([qpd_data.YLatt_ave],15)
    xlabel('Y Latt. Ave. (V)')
    
    axLYSD=subplot(3,4,8);
    histogram([qpd_data.YLatt_std],15)
    xlabel('Y Latt. Std. Dev. (V)')
    
    axLZA=subplot(3,4,9);
    histogram([qpd_data.ZLatt_ave],15)
    xlabel('Z Latt. Ave. (V)')
    
    axLZSD=subplot(3,4,10);
    histogram([qpd_data.ZLatt_std],15)
    xlabel('Z Latt. Std. Dev. (V)')

end

%% Save the figures and QPD analysis

if opts.doSave

    saveas(fmsum,[saveDir filesep 'qpd_modulation_summary.png']);
    saveas(fmsum_um,[saveDir filesep 'qpd_modulation_summary_um.png']);
    saveas(fpsum,[saveDir filesep 'qpd_power_summary.png']);

    filename=fullfile(saveDir,'qpddata.mat');
    save(filename,'qpd_out');

    
end


end