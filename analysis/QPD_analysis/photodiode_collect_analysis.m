function [pd_summary,f_trace,f_mod,f_sum] = photodiode_collect_analysis(pd_data,opts)

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analyzes list of QPD filenames 
% Outputs average powers, fit parameters, raw channels with time shifted to
% start of modulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Collecting total analysis for QPD files');

%% Assign analyzed data to output

pd_summary.pd_data = pd_data;

%% Compute average parameters for all runs
% Take average and standard deviation of fit values during first 50ms
    % of modulation
pd_summary.mod_params.X1.A     = mean([pd_data.modfit_X1_A]);
pd_summary.mod_params.X1.A_err = std([pd_data.modfit_X1_A]);
pd_summary.mod_params.X1.B     = mean([pd_data.modfit_X1_B]);
pd_summary.mod_params.X1.B_err = std([pd_data.modfit_X1_B]);
pd_summary.mod_params.X1.C     = mean([pd_data.modfit_X1_C]);
pd_summary.mod_params.X1.C_err = std([pd_data.modfit_X1_C]);
pd_summary.mod_params.X1.T     = mean([pd_data.modfit_X1_T]);
pd_summary.mod_params.X1.T_err = std([pd_data.modfit_X1_T]);

pd_summary.mod_params.X2.A     = mean([pd_data.modfit_X2_A]);
pd_summary.mod_params.X2.A_err = std([pd_data.modfit_X2_A]);
pd_summary.mod_params.X2.B     = mean([pd_data.modfit_X2_B]);
pd_summary.mod_params.X2.B_err = std([pd_data.modfit_X2_B]);
pd_summary.mod_params.X2.C     = mean([pd_data.modfit_X2_C]);
pd_summary.mod_params.X2.C_err = std([pd_data.modfit_X2_C]);
pd_summary.mod_params.X2.T     = mean([pd_data.modfit_X2_T]);
pd_summary.mod_params.X2.T_err = std([pd_data.modfit_X2_T]);

% Take average and standard deviation of average powers during -350ms to 50ms
pd_summary.powers.ODT1      = mean([pd_data.ODT1_ave]);
pd_summary.powers.ODT1_err  = std([pd_data.ODT1_ave]);
pd_summary.powers.ODT2      = mean([pd_data.ODT2_ave]);
pd_summary.powers.ODT2_err  = std([pd_data.ODT2_ave]);
pd_summary.powers.XLatt     = mean([pd_data.XLatt_ave]);
pd_summary.powers.XLatt_err = std([pd_data.XLatt_ave]);
pd_summary.powers.YLatt     = mean([pd_data.YLatt_ave]);
pd_summary.powers.YLatt_err = std([pd_data.YLatt_ave]);
pd_summary.powers.ZLatt     = mean([pd_data.ZLatt_ave]);
pd_summary.powers.ZLatt_err = std([pd_data.ZLatt_ave]);

%% 


f_trace = figure;
f_trace.Position = [100 100 540 400]; 
clf
f_trace.Color= 'w';

cmap =jet(length(pd_data));

% Plot all data (Not sure what the purpose of this figure is?
for nn=1:length(pd_data)
tt=linspace(min(pd_data(nn).modfit_t),max(pd_data(nn).modfit_t),1e3);
%       % Plot X1
        subplot(211)
        plot(pd_data(nn).modfit_t,pd_data(nn).modfit_normX1data,'-','color',cmap(nn,:))
        hold on;
        plot(tt,feval(pd_data(nn).modfit_X1,tt),'-','color',cmap(nn,:))
        xlim([0 50])
        xlabel('Time (ms)')
        ylabel('Normalized X1')

        % Plot X2
        subplot(212)
        plot(pd_data(nn).modfit_t,pd_data(nn).modfit_normX2data,'-','color',cmap(nn,:))
        hold on;
        plot(tt,feval(pd_data(nn).modfit_X2,tt),'color',cmap(nn,:))    
        xlim([0 50])
        xlabel('Time (ms)')
        ylabel('Normalized X2')        
end

%% Plot All Traces and histgorams

f_mod = figure;
f_mod.Position = [100 100 1200 800]; 
clf
f_mod.Color= 'w';
plotColor = jet(length(pd_data));


% Select the data during the modulation, from ramp up to end
% ramp_start = -200;
% ramp_start = find(round(pd_data(nn).t,1)==-ixondata(nn).Params.conductivity_mod_ramp_time);

xL = [-200 150];
for nn=1:length(pd_data)
        %%%%%%%%%%%% Plot summary of modulation traces %%%%%%%%%%%%%%

        % Find end of mod from lattice ramp up BUGGY RL FIX
        % (approx. as max slope of lattice ramp)
%         mod_end   = find(diff(pd_data(nn).XLATT)==max(diff(pd_data(nn).XLATT)));
%         mod_end   = max(find(diff(pd_data(nn).TRIG)==-1));

        % Add 50 ms buffer time to each end 
%         buffer_t = 250;         

        ax1=subplot(4,3,1);
        plot(pd_data(nn).t,pd_data(nn).X1./pd_data(nn).SUM1,'.',  'color',  plotColor(nn,:))
        hold on;
        xlabel('Time (ms)')
        ylabel('ODT 1 X (V/V)')
    
        ax2=subplot(4,3,4);
        plot(pd_data(nn).t,pd_data(nn).Y1./pd_data(nn).SUM1,'.',  'color',  plotColor(nn,:))
        hold on;
        xlabel('Time (ms)')
        ylabel('ODT 1 Y (V/V)')
    
        ax3=subplot(4,3,7);
        plot(pd_data(nn).t,pd_data(nn).X2./pd_data(nn).SUM2,'.',  'color',  plotColor(nn,:))
        hold on;
        xlabel('Time (ms)')
        ylabel('ODT 2 X (V/V)')
    
        ax4=subplot(4,3,10);
        plot(pd_data(nn).t,pd_data(nn).Y2./pd_data(nn).SUM2,'.',  'color',  plotColor(nn,:))
        hold on;
        xlabel('Time (ms)')
        ylabel('ODT 2 Y (V/V)')
    
        linkaxes([ax1 ax2 ax3 ax4],'x')
%         xlim([min(pd_data(nn).t(ramp_start-buffer_t:mod_end+buffer_t)) max(pd_data(nn).t(ramp_start-buffer_t:mod_end+buffer_t))])
        xlim(xL);
  ylim('auto')    
end

    ax5=subplot(4,3,[2]);
    pd = fitdist([pd_data.modfit_X1_A]','Normal');
    a=histfit([pd_data.modfit_X1_A],15);
    str = ['$' num2str(pd.mu,'%.2e') ' \pm ' num2str(pd.sigma,'%.2e') '$'];
    legend(a(2),{str},'interpreter','latex');
    xlabel('X1/SUM1 amplitude (V/V)')

    ax5b=subplot(4,3,[5]);
    pd = fitdist([pd_data.modfit_X1_C]','Normal');
    a=histfit([pd_data.modfit_X1_C],15);
    str = ['$' num2str(pd.mu,'%.2e') ' \pm ' num2str(pd.sigma,'%.2e') '$'];
    legend(a(2),{str},'interpreter','latex');
    xlabel('X1/SUM1 offset (V/V)')
    
    ax6=subplot(4,3,[8]);
    pd = fitdist([pd_data.modfit_X2_A]','Normal');
    a=histfit([pd_data.modfit_X2_A],15);
    str = ['$' num2str(pd.mu,'%.2e') ' \pm ' num2str(pd.sigma,'%.2e') '$'];
    legend(a(2),{str},'interpreter','latex');
    xlabel('X2/SUM2 amplitude (V/V)')
    
    ax6b=subplot(4,3,[11]);
    pd = fitdist([pd_data.modfit_X2_C]','Normal');
    a=histfit([pd_data.modfit_X2_C],15);
    str = ['$' num2str(pd.mu,'%.2e') ' \pm ' num2str(pd.sigma,'%.2e') '$'];
    legend(a(2),{str},'interpreter','latex');
    xlabel('X2/SUM2 offset (V/V)')

    ax7=subplot(4,3,[3 6]);
    histogram([pd_data.modfit_X1_B]/(2*pi),15)
    hold on;
    histogram([pd_data.modfit_X2_B]/(2*pi),15)
    xlabel('Phase/ (2$\pi$)','Interpreter','Latex','FontName','Helvetica')
legend({'odt1','odt2'});
    
    ax8=subplot(4,3,[9 12]);
    histogram(1./[pd_data.modfit_X1_T],15)
    hold on;
    histogram(1./[pd_data.modfit_X2_T],15)
    xlabel('Frequency (kHz)')
legend({'odt1','odt2'});

    ax7=subplot(4,3,[3 6]);



%%
%     figure(201)

% How necessaryis this?


        
        % Plot modulation in um

%         fmsum_um=figure(401);
%         fmsum_um.Position = [600 100 800 600];
%     
%         ax1=subplot(4,3,1);
%         plot(pd_data(nn).t,pd_data(nn).X1./pd_data(nn).SUM1*94.0311,'.',  'color',  plotColor(nn,:))
%         hold on;
%         xlabel('Time (ms)')
%         ylabel('ODT 1 X (um)')
%     
%         ax2=subplot(4,3,4);
%         plot(pd_data(nn).t,pd_data(nn).Y1./pd_data(nn).SUM1,'.',  'color',  plotColor(nn,:))
%         hold on;
%         xlabel('Time (ms)')
%         ylabel('ODT 1 Y (V/V)')
%     
%         ax3=subplot(4,3,7);
%         plot(pd_data(nn).t,pd_data(nn).X2./pd_data(nn).SUM2*994.9276,'.',  'color',  plotColor(nn,:))
%         hold on;
%         xlabel('Time (ms)')
%         ylabel('ODT 2 X (um)')
%     
%         ax4=subplot(4,3,10);
%         plot(pd_data(nn).t,pd_data(nn).Y2./pd_data(nn).SUM2,'.',  'color',  plotColor(nn,:))
%         hold on;
%         xlabel('Time (ms)')
%         ylabel('ODT 2 Y (V/V)')
%     
%         linkaxes([ax1 ax2 ax3 ax4],'x')
% %         xlim([min(pd_data(nn).t(ramp_start-buffer_t:mod_end+buffer_t)) max(pd_data(nn).t(ramp_start-buffer_t:mod_end+buffer_t))])
%         ylim('auto')    
   
% 
%% Summary plot for the modulation %%%%%%%%%%%%%%%
% 
%     % um
% How necessaryis this?
%     figure(401)
% 
%     % um Feb 13/24 conversion 94.0311 um/(V/V)
%     ax5=subplot(4,3,[2 5]);
%     histogram([pd_data.modfit_X1_A]*94.0311,15)
%     xlabel('Norm. X1 mod. amp. (um)')
%     
%     % um Feb 13/24 conversion 994.9276 um/(V/V)
%     ax6=subplot(4,3,[3 6]);
%     histogram([pd_data.modfit_X2_A]*994.9276,15)
%     xlabel('Norm. X2 mod. amp. (um)')
%     
%     ax7=subplot(4,3,[8 11]);
%     histogram([pd_data.modfit_X1_B]/(2*pi),15)
%     hold on;
%     histogram([pd_data.modfit_X2_B]/(2*pi),15)
%     xlabel('Phase/ (2$\pi$)','Interpreter','Latex','FontName','Helvetica')
%     
%     ax8=subplot(4,3,[9 12]);
%     histogram(1./[pd_data.modfit_X1_T],15)
%     hold on;
%     histogram(1./[pd_data.modfit_X2_T],15)
%     xlabel('Frequency (kHz)')
    
%%
%     %%%%%%%%%%%%% Summary plot for the power %%%%%%%%%%%%%
%


    f_sum = figure;
    f_sum.Position = [100 100 1200 800]; 
    clf
    f_sum.Color= 'w';
    axODT1A=subplot(3,4,1);
    histogram([pd_data.ODT1_ave],15)
    xlabel('Sum 1 Ave. (V)')
    
    axODT1SD=subplot(3,4,2);
    histogram([pd_data.ODT1_std],15)
    xlabel('Sum 1 Std. Dev. (V)')
    
    axODT2A=subplot(3,4,3);
    histogram([pd_data.ODT2_ave],15)
    xlabel('Sum 2 Ave. (V)')
    
    axODT2SD=subplot(3,4,4);
    histogram([pd_data.ODT2_std],15)
    xlabel('Sum 2 Std. Dev. (V)')
    
    axLXA=subplot(3,4,5);
    histogram([pd_data.XLatt_ave],15)
    xlabel('X Latt. Ave. (V)')
    
    axLXSD=subplot(3,4,6);
    histogram([pd_data.XLatt_std],15)
    xlabel('X Latt. Std. Dev. (V)')
    
    axLYA=subplot(3,4,7);
    histogram([pd_data.YLatt_ave],15)
    xlabel('Y Latt. Ave. (V)')
    
    axLYSD=subplot(3,4,8);
    histogram([pd_data.YLatt_std],15)
    xlabel('Y Latt. Std. Dev. (V)')
    
    axLZA=subplot(3,4,9);
    histogram([pd_data.ZLatt_ave],15)
    xlabel('Z Latt. Ave. (V)')
    
    axLZSD=subplot(3,4,10);
    histogram([pd_data.ZLatt_std],15)
    xlabel('Z Latt. Std. Dev. (V)')

 


end