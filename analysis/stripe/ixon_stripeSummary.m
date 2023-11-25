function hF2 = ixon_stripeSummary(stripes,xvals,opts)

global ixon_imgdir


%% Process Fits
[cface1,cedge1] = ixoncolororder(1);

% Summarize the results
hF2=figure;
hF2.Position=[100 100 1000 500];
hF2.Color='w';
clf;

uicontrol('style','text','string','iXon, stripe','units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',12,'fontweight','bold',...
    'position',[2 2 100 20]);

% Folder directory
strs=strsplit(ixon_imgdir,filesep);
str=[strs{end-1} filesep strs{end}];
% Image directory folder string
t=uicontrol('style','text','string',str,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',6);
t.Position(4)=t.Extent(4);
t.Position(3)=hF2.Position(3);
t.Position(1:2)=[5 hF2.Position(4)-t.Position(4)];


% Wavelength
subplot(221);
errorbar(xvals,[stripes.L],[stripes.L_err],'marker','o',...
    'MarkerFacecolor',cface1,'markeredgecolor',cedge1,'linestyle','none',...
    'linewidth',1.5,'color',cedge1);
xlabel([opts.xVar ' (' opts.xUnit ')'],'interpreter','none');
ylabel('wavelength (px)');
grid on

if isequal(opts.xVar,'ExecutionDate')
    datetick('x');
    xlabel('time','fontsize',10);
end
xlim([min(xvals) max(xvals)]);

% Angle
axb2=subplot(222);
errorbar(xvals,[stripes.theta],[stripes.theta_err],'marker','o',...
    'MarkerFacecolor',cface1,'markeredgecolor',cedge1,'linestyle','none',...
    'linewidth',1.5,'color',cedge1);
xlabel([opts.xVar ' (' opts.xUnit ')'],'interpreter','none');
ylabel('angle (deg.)');
grid on
if isequal(opts.xVar,'ExecutionDate')
    datetick('x');
    xlabel('time','fontsize',10);
end
xlim([min(xvals) max(xvals)]);

% Phase
subplot(223);
errorbar(xvals,[stripes.phi]/(2*pi),[stripes.phi_err]/(2*pi),'marker','o',...
    'MarkerFacecolor',cface1,'markeredgecolor',cedge1,'linestyle','none',...
    'linewidth',1.5,'color',cedge1);
xlabel([opts.xVar ' (' opts.xUnit ')'],'interpreter','none');
ylabel('unwrapped phase (2\pi)');
grid on
if isequal(opts.xVar,'ExecutionDate')
    datetick('x');
    xlabel('time','fontsize',10);
end
xlim([min(xvals) max(xvals)]);

% mod depth
subplot(224);
errorbar(xvals,[stripes.B],[stripes.B_err],'marker','o',...
    'MarkerFacecolor',cface1,'markeredgecolor',cedge1,'linestyle','none',...
    'linewidth',1.5,'color',cedge1);
xlabel([opts.xVar ' (' opts.xUnit ')'],'interpreter','none');
ylabel('modulatoin depth');
grid on
if isequal(opts.xVar,'ExecutionDate')
    datetick('x');
    xlabel('time','fontsize',10);
end
xlim([min(xvals) max(xvals)]);
end

