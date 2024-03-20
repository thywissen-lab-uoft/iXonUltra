function [pd_data,lattice_summary,hF]=pd_lattice_load(pd_data,opts)

hF=[];

hF = figure;
hF.Color = 'w';
hF.Name = 'photiode lattice load';
hF.Position = [100 100 750 500];
%% Analyze the lattice load traces
c=jet(length(pd_data));
for nn = 1:length(pd_data)    
    % Period of 60Hz in ms
    T60 = 1e3/60;
    
    t = pd_data(nn).t;
    v = pd_data(nn).v;
    
    % Indeces before modulation
    ii = t<=0;

    % Get the lattice photodiode information
    xlatt = v(:,7);
    ylatt = v(:,8);
    zlatt = v(:,9);

    % Get the photodiode information before loading
    tpre  = t(ii);
    xlatt = xlatt(ii);
    ylatt = ylatt(ii);
    zlatt = zlatt(ii);
    
    % How many 60 Hz periods to measure the zero level
    tzero = T60*6;    
    
    iL = find(tpre<=(tzero+tpre(1)));
    
    % Measure the zero level mean
    xlatt_L_mean = mean(xlatt(iL));
    ylatt_L_mean = mean(ylatt(iL));
    zlatt_L_mean = mean(zlatt(iL));

    % Measure the zero level standard deviation
    xlatt_L_std = std(xlatt(iL));
    ylatt_L_std = std(ylatt(iL));
    zlatt_L_std = std(zlatt(iL));
    
    % How many 60 Hz periods to measure the high level
    thigh = T60*6;    
    
    iH = find([tpre<=0].*[tpre>=-thigh]);
    
    % Measure the zero level mean
    xlatt_H_mean = mean(xlatt(iH));
    ylatt_H_mean = mean(ylatt(iH));
    zlatt_H_mean = mean(zlatt(iH));

    % Measure the zero level standard deviation
    xlatt_H_std = std(xlatt(iH));
    ylatt_H_std = std(ylatt(iH));
    zlatt_H_std = std(zlatt(iH));
    
    % Measure the x lattice delta
    xlatt_val = xlatt_H_mean-xlatt_L_mean;
    xlatt_err = sqrt(xlatt_H_std.^2+xlatt_L_std.^2);

    % Measure the y lattice delta
    ylatt_val = ylatt_H_mean-ylatt_L_mean;
    ylatt_err = sqrt(ylatt_H_std.^2+ylatt_L_std.^2);
    
    % Measure the z lattice delta
    zlatt_val = zlatt_H_mean-zlatt_L_mean;
    zlatt_err = sqrt(zlatt_H_std.^2+zlatt_L_std.^2);
    
    pd_data(nn).xlatt_val = xlatt_val;
    pd_data(nn).xlatt_err = xlatt_err;
    pd_data(nn).ylatt_val = ylatt_val;
    pd_data(nn).ylatt_err = ylatt_err;
    pd_data(nn).zlatt_val = zlatt_val;
    pd_data(nn).zlatt_err = zlatt_err;    
    
    figure(hF);
    ax1=subplot(3,3,[1 2]);
    plot(tpre,xlatt-xlatt_L_mean,'-','linewidth',.5);
    hold on
    ylabel('xlatt - xlatt0 (mV)');   
    title('x lattice');
    xlim([tpre(1) tpre(end)]);
    
    ax2=subplot(3,3,[4 5]);
    plot(tpre,ylatt-ylatt_L_mean,'-','linewidth',.5);
    hold on
    ylabel('ylatt - ylatt0 (mV)');   
    title('y lattice');
    xlim([tpre(1) tpre(end)]);

    ax3=subplot(3,3,[7 8]);
    plot(tpre,zlatt-zlatt_L_mean,'-','linewidth',.5);
    hold on
    ylabel('zlatt - zlatt0 (mV)');   
    title('z lattice');
    xlim([tpre(1) tpre(end)]);
    xlabel('times (ms)');
end

if nargin >1 && isfield(opts,'FigLabel') && ~isempty(opts.FigLabel)        
    tFig=uicontrol('style','text','string',opts.FigLabel,...
        'units','pixels','backgroundcolor',...
        'w','horizontalalignment','left');
    tFig.Position(4)=tFig.Extent(4);
    tFig.Position(3)=hF.Position(3);
    tFig.Position(1:2)=[5 hF.Position(4)-tFig.Position(4)];
end

%% Plot
    % X Voltages
    vX = mean([pd_data.xlatt_val]);
    eX = mean([pd_data.xlatt_err]);
    strX = [num2str(round(vX,1)) '\pm' num2str(round(eX,1)) ' mV'];    
    % X Calibration
    cX = opts.XCalibration;
    cXerr = opts.XCalibrationUncertainty;    
    % X Depth
    uX = vX*cX;
    uXerr = uX*sqrt((eX/vX).^2+(cXerr/cX).^2);
    
    
    % Y Voltages
    vY = mean([pd_data.ylatt_val]);
    eY = mean([pd_data.ylatt_err]);
    strY = [num2str(round(vY,1)) '\pm' num2str(round(eY,1)) ' mV'];
    % Y Calibration
    cY = opts.YCalibration;
    cYerr = opts.YCalibrationUncertainty;    
    % Y Depth
    uY = vY*cY;
    uYerr = uY*sqrt((eY/vY).^2+(cYerr/cY).^2);
    
    % Z Voltages
    vZ = mean([pd_data.zlatt_val]);
    eZ = mean([pd_data.zlatt_err]);
    strZ = [num2str(round(vZ,1)) '\pm' num2str(round(eZ,1)) ' mV'];
    % Z Calibration
    cZ = opts.ZCalibration;
    cZerr = opts.ZCalibrationUncertainty;    
    % Z Depth
    uZ = vZ*cZ;
    uZerr = uZ*sqrt((eZ/vZ).^2+(cZerr/cZ).^2);
    
    
    %%
    
    axes(ax1)
    plot([min(tpre) max(tpre)],[1 1]*vX,'-','color',[.3 .3 .3]);
    plot([min(tpre) max(tpre)],[0 0],'-','color',[.3 .3 .3]);   
    arrow('Start',[min(tpre)+20 vX/2],'Stop',[min(tpre)+20 vX],'length',10);
    arrow('Start',[min(tpre)+20 vX/2],'Stop',[min(tpre)+20 0],'length',10);    
    text(min(tpre)+40,mean([pd_data.xlatt_val]),strX,'horizontalalignment','left',...
        'verticalalignment','top');    
    
    axes(ax2)
    plot([min(tpre) max(tpre)],[1 1]*vY,'-','color',[.3 .3 .3]);
    plot([min(tpre) max(tpre)],[0 0],'-','color',[.3 .3 .3]);   
    arrow('Start',[min(tpre)+20 vY/2],'Stop',[min(tpre)+20 vY],'length',10);
    arrow('Start',[min(tpre)+20 vY/2],'Stop',[min(tpre)+20 0],'length',10);    
    text(min(tpre)+40,mean([pd_data.ylatt_val]),strY,'horizontalalignment','left',...
        'verticalalignment','top');   
    axes(ax3)
    plot([min(tpre) max(tpre)],[1 1]*vZ,'-','color',[.3 .3 .3]);
    plot([min(tpre) max(tpre)],[0 0],'-','color',[.3 .3 .3]);   
    arrow('Start',[min(tpre)+20 vZ/2],'Stop',[min(tpre)+20 vZ],'length',10);
    arrow('Start',[min(tpre)+20 vZ/2],'Stop',[min(tpre)+20 0],'length',10);    
    text(min(tpre)+40,mean([pd_data.zlatt_val]),strZ,'horizontalalignment','left',...
        'verticalalignment','top');
    
%% Tables
    % char(945:969)
    % char(913:937)
    
    axt1=subplot(3,3,3);
    t1 = uitable('units','normalized','RowName',{},'ColumnName',{'name','value','error'},...
        'ColumnWidth',{75 65 65});
    t1.Position(1:2) = axt1.Position(1:2);
    delete(axt1);    
    d1 = {[char(916) 'xlatt (mV)'], num2str(round(vX,1)), num2str(round(eX,1));
        'calib (Er/mV)', num2str(cX), num2str(cXerr);
        'calib date', opts.XCalibrationString,'';
        'xlatt (Er)', num2str(round(uX,2)), num2str(round(uXerr,2))};    
    t1.Data=d1;
    t1.Position(3:4)=t1.Extent(3:4);

    
    axt2=subplot(3,3,6);
    t2 = uitable('units','normalized','RowName',{},'ColumnName',{'name','value','error'},...
        'ColumnWidth',{75 65 65});
    t2.Position(1:2) = axt2.Position(1:2);
    delete(axt2);    
    d2 = {[char(916) 'ylatt (mV)'], num2str(round(vY,1)), num2str(round(eY,1));
        'calib (Er/mV)', num2str(cY), num2str(cYerr);
        'calib date', opts.YCalibrationString,'';
        'ylatt (Er)', num2str(round(uY,2)), num2str(round(uYerr,2))};    
    t2.Data=d2;
    t2.Position(3:4)=t2.Extent(3:4);
    
        
    axt3=subplot(3,3,9);
    t3 = uitable('units','normalized','RowName',{},'ColumnName',{'name','value','error'},...
        'ColumnWidth',{75 65 65});
    t3.Position(1:2) = axt3.Position(1:2);
    delete(axt3);    
    d3 = {[char(916) 'zlatt (mV)'], num2str(round(vZ,1)), num2str(round(eZ,1));
        'calib (Er/mV)', num2str(cZ), num2str(cZerr);
        'calib date', opts.ZCalibrationString,'';
        'zlatt (Er)', num2str(round(uZ,2)), num2str(round(uZerr,2))};    
    t3.Data=d3;
    t3.Position(3:4)=t3.Extent(3:4);
    
    %%
    
    lattice_summary = struct;    
    lattice_summary.xlatt_Er = uX;
    lattice_summary.xlatt_Er_err = uXerr;
    lattice_summary.ylatt_Er = uY;
    lattice_summary.ylatt_Er_err = uYerr;  
    lattice_summary.zlatt_Er = uZ;
    lattice_summary.zlatt_Er_err = uZerr;
    
    lattice_summary.xlatt_mV = vX;
    lattice_summary.xlatt_mV_err = eX;
    lattice_summary.xlatt_calib= cX;
    lattice_summary.xlatt_calib_err= cXerr;
    lattice_summary.xlatt_calib_str= opts.XCalibrationString;

    lattice_summary.ylatt_mV = vY;
    lattice_summary.ylatt_mV_err = eY;
    lattice_summary.ylatt_calib= cY;
    lattice_summary.ylatt_calib_err= cYerr;
    lattice_summary.ylatt_calib_str= opts.YCalibrationString;
    
    lattice_summary.zlatt_mV = vZ;
    lattice_summary.zlatt_mV_err = eZ;
    lattice_summary.zlatt_calib= cZ;
    lattice_summary.zlatt_calib_err= cZerr;
    lattice_summary.zlatt_calib_str= opts.ZCalibrationString;
end

