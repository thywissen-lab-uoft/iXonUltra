function [pd_data,modulation_summary,hF] = pd_odt_modulation(pd_data,opts)

hF=[];

hF = figure;
hF.Color = 'w';
hF.Name = 'photiode odt modulation';
hF.Position = [100 100 1400 500];
%% Analyze the lattice load traces

% odt1all = [];
% odt2all = [];
c=jet(length(pd_data));
for nn = 1:length(pd_data)      
    t = pd_data(nn).t;
    v = pd_data(nn).v;      
    
    % Modulation frequency in kHz
    freq = 1e-3*pd_data(nn).Params.conductivity_mod_freq;    
    % Ramp time
    t0  = pd_data(nn).Params.conductivity_mod_ramp_time; 
    % Modulation time
    tmod = pd_data(nn).Params.conductivity_mod_time;

    % Get the odt photodiode information
    x1 = v(:,1)./v(:,3);
    x2 = v(:,4)./v(:,6);

    % Time limits
    tlim = [0 tmod];
    ia = find(t>=(tlim(1)),1);
    ib = find(t>=(tlim(2)),1);
    
    % Truncated data
    tsub = t(ia:ib)+t0;
    x1sub = x1(ia:ib);
    x2sub = x2(ia:ib);        
    
     % Basic sin func to fit
    myfunc = @(xA,x0,phi,t) -xA*sin(2*pi*freq*t+phi) + x0;    
    myfit = fittype(@(xA,x0,phi,t) myfunc(xA,x0,phi,t),'independent',{'t'},...
        'coefficients',{'xA','x0','phi'});
    opt = fitoptions(myfit);

    % ODT1 guess
    xAg1 = (max(x1sub)-min(x1sub))/2;
    x0g1 = (max(x1sub)+min(x1sub))/2;

    % ODT2 guess
    xAg2 = (max(x2sub)-min(x2sub))/2;
    x0g2 = (max(x2sub)+min(x2sub))/2;
    
    opt.StartPoint = [xAg1 x0g1 0];    
    fout1 = fit(tsub',x1sub,myfit,opt);  
  
    opt.StartPoint = [xAg2 x0g2 0];    
    fout2 = fit(tsub',x2sub,myfit,opt);     

    figure(hF);
    ax1=subplot(2,6,[1 2]);
    plot(t+t0,x1,'.','linewidth',.5,'color',c(nn,:));
    hold on
    ylabel('odt1 x (mV/mV)');   
    title('odt1');
    xlim([0 t0+tmod]);    
    plot(tsub,feval(fout1,tsub),'-','color',c(nn,:));    
    xlabel('total modulation time (ms)');
    
    ax2=subplot(2,6,[7 8]);
    plot(t+t0,x2,'.','linewidth',.5,'color',c(nn,:));
    hold on
    ylabel('odt2 x (mV/mV)');   
    title('odt2');    
    xlim([0 t0+tmod]);    
    xlabel('time (ms)');
    plot(tsub,feval(fout2,tsub),'-','color',c(nn,:));
        xlabel('total modulation time (ms)');

    c1 = confint(fout1,0.667);
    c2 = confint(fout2,0.667);

    pd_data(nn).x1_amplitude = fout1.xA;
    pd_data(nn).x1_amplitude_err = 0.5*(c1(2,1)-c1(1,1));
    pd_data(nn).x1_offset = fout1.x0;
    pd_data(nn).x1_offset_err = 0.5*(c1(2,2)-c1(1,2));
    pd_data(nn).x1_phase = fout1.phi;
    pd_data(nn).x1_phase_err = 0.5*(c1(2,3)-c1(1,3));
    
    
    pd_data(nn).x2_amplitude = fout2.xA;
    pd_data(nn).x2_amplitude_err = 0.5*(c2(2,1)-c2(1,1));
    pd_data(nn).x2_offset = fout2.x0;
    pd_data(nn).x2_offset_err = 0.5*(c2(2,2)-c2(1,2));    
    pd_data(nn).x2_phase = fout2.phi;
    pd_data(nn).x2_phase_err = 0.5*(c2(2,3)-c2(1,3));

    
    pd_data(nn).x1fit = fout1;
    pd_data(nn).x2fit = fout2;
end

p       = [pd_data.Params];
tmod    = [p.conductivity_mod_time];
co      = get(gca,'colororder');

subplot(2,6,3,'parent',hF);
errorbar(tmod,[pd_data.x1_amp],[pd_data.x1_amp_err],'o','color',.5*co(1,:),...
    'markerfacecolor',co(1,:),'linewidth',1,'markeredgecolor',co(1,:)*.5);
xlabel('conductivity mod time (ms)');
ylabel('amplitude (mV/mV)');
title('odt1 amplitude');

subplot(2,6,4,'parent',hF);
errorbar(tmod,[pd_data.x1_offset],[pd_data.x1_offset_err],'o','color',.5*co(1,:),...
    'markerfacecolor',co(1,:),'linewidth',1,'markeredgecolor',co(1,:)*.5);
xlabel('conductivity mod time (ms)');
ylabel('offset (mV/mV)');
title('odt1 offset');


subplot(2,6,5,'parent',hF);
errorbar(tmod,[pd_data.x1_phase],[pd_data.x1_phase_err],'o','color',.5*co(1,:),...
    'markerfacecolor',co(1,:),'linewidth',1,'markeredgecolor',co(1,:)*.5);
xlabel('conductivity mod time (ms)');
ylabel('phase (rad)');
title('odt1 phase');



subplot(2,6,9,'parent',hF);
errorbar(tmod,[pd_data.x2_amp],[pd_data.x2_amp_err],'o','color',.5*co(2,:),...
    'markerfacecolor',co(2,:),'linewidth',1,'markeredgecolor',co(2,:)*.5);
xlabel('conductivity mod time (ms)');
ylabel('amplitude (mV/mV)');
title('odt2 amplitude');

subplot(2,6,10,'parent',hF);
errorbar(tmod,[pd_data.x2_offset],[pd_data.x2_offset_err],'o','color',.5*co(2,:),...
    'markerfacecolor',co(2,:),'linewidth',1,'markeredgecolor',co(2,:)*.5);
xlabel('conductivity mod time (ms)');
ylabel('offset (mV/mV)');
title('odt2 offset');

subplot(2,6,11,'parent',hF);
errorbar(tmod,[pd_data.x2_phase],[pd_data.x2_phase_err],'o','color',.5*co(2,:),...
    'markerfacecolor',co(2,:),'linewidth',1,'markeredgecolor',co(2,:)*.5);
xlabel('conductivity mod time (ms)');
ylabel('phase (rad)');
title('odt2 phase');

%% Process

x1A = mean([pd_data.x1_amplitude]);
x1A_err = std([pd_data.x1_amplitude]);
x10 = mean([pd_data.x1_offset]);
x10_err = std([pd_data.x1_offset]);
x1phi = mean([pd_data.x1_phase]);
x1phi_err = std([pd_data.x1_phase]);

x2A = mean([pd_data.x2_amplitude]);
x2A_err = std([pd_data.x2_amplitude]);
x20 = mean([pd_data.x2_offset]);
x20_err = std([pd_data.x2_offset]);
x2phi = mean([pd_data.x2_phase]);
x2phi_err = std([pd_data.x2_phase]);

modulation_summary = struct;

dR = opts.ODT1Calibration*x1A + ...
    opts.ODT2Calibration*x2A;

R0 = opts.ODT1Calibration*x10 + ...
    opts.ODT2Calibration*x20;

% Assume 5% calibration error (arbitrary)
dR_err = sqrt((x1A_err*opts.ODT1Calibration).^2 + ...
    (x1A*opts.ODT1Calibration*0.05).^2 + ...
    (x2A_err*opts.ODT2Calibration).^2 + ...
    (x2A*opts.ODT2Calibration*0.05).^2);

R0_err = sqrt((x10_err*opts.ODT1Calibration).^2 + ...
    (x10*opts.ODT1Calibration*0.05).^2 + ...
    (x20_err*opts.ODT2Calibration).^2 + ...
    (x20*opts.ODT2Calibration*0.05).^2);

modulation_summary.xlattice_amplitude_um = dR(1);
modulation_summary.xlattice_amplitude_um_err = dR_err(1);

modulation_summary.ylattice_amplitude_um = dR(2);
modulation_summary.ylattice_amplitude_um_err = dR_err(2);

modulation_summary.xlattice_offset_um = R0(1);
modulation_summary.xlattice_offset_um_err = R0_err(1);

modulation_summary.ylattice_offset_um = R0(2);
modulation_summary.ylattice_offset_um_err = R0_err(2);

modulation_summary.phase_rad = mean([x1phi x2phi]);
modulation_summary.freq_Hz =pd_data(1).Params.conductivity_mod_freq;

modulation_summary.x1_amplitude = x1A;
modulation_summary.x1_amplitude_err = x1A_err;
modulation_summary.x1_offset = x10;
modulation_summary.x1_offset_err = x10_err;
modulation_summary.x1_phase = x1phi;
modulation_summary.x1_phase_err = x1phi_err;


modulation_summary.x2_amplitude = x2A;
modulation_summary.x2_amplitude_err = x2A_err;
modulation_summary.x2_offset = x20;
modulation_summary.x2_offset_err = x20_err;
modulation_summary.x2_phase = x2phi;
modulation_summary.x2_phase_err = x2phi_err;

%%



axt = subplot(2,6,[6 12],'parent',hF);
myt = uitable('units','normalized','RowName',{},'ColumnName',{'name','value','error'},...
    'ColumnWidth',{95 65 65});
myt.Position(1:2) = axt.Position(1:2);
delete(axt);    
d1 = {'xlatt amp (um)', num2str(round(dR(1),1)), num2str(round(dR_err(1),2));
    'ylatt amp (um)', num2str(round(dR(2),1)), num2str(round(dR_err(2),2));
    'phase (rad)',  num2str(mean([x1phi x2phi])),'';
    'xlatt offset (um)', num2str(round(R0(1),1)), num2str(round(R0_err(1),1));
    'ylatt offset (um)', num2str(round(R0(2),1)), num2str(round(R0_err(2),1));    
    'odt calib', opts.ODTCalibrationString,'';
    'odt1 x (um/V/V)', num2str(round(opts.ODT1Calibration(1))), '';
    'odt1 y (um/V/V)', num2str(round(opts.ODT1Calibration(2))), '';
    'odt2 x (um/V/V)', num2str(round(opts.ODT2Calibration(1))), '';
    'odt2 y (um/V/V)', num2str(round(opts.ODT2Calibration(2))), ''};
myt.Data=d1;
myt.Position(3:4)=myt.Extent(3:4);
    

if nargin >1 && isfield(opts,'FigLabel') && ~isempty(opts.FigLabel)        
    tFig=uicontrol('style','text','string',opts.FigLabel,...
        'units','pixels','backgroundcolor',...
        'w','horizontalalignment','left','parent',hF);
    tFig.Position(4)=tFig.Extent(4);
    tFig.Position(3)=hF.Position(3);
    tFig.Position(1:2)=[5 hF.Position(4)-tFig.Position(4)];
end



end

