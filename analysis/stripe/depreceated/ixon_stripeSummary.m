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


% Phase
fitWavelength=1;
if fitWavelength
    
    myfit = fittype('L/(I-I0)',...
        'coefficients',{'L','I0'},...
        'independent','I');

    % Fit options and guess
    opt=fitoptions(myfit);        
    Lg = 300;
    I0g = -1;
    G=[Lg I0g];        
    opt.StartPoint=G;

    % Perform the fit
    fout=fit(xvals',[stripes.L]',myfit,opt)

    l_text = sprintf('\\lambda =%.2f px/(I - (%.2f A))',fout.L,fout.I0);
    end

% Wavelength
subplot(221);
errorbar(xvals,[stripes.L],[stripes.L_err],'marker','o',...
    'MarkerFacecolor',cface1,'markeredgecolor',cedge1,'linestyle','none',...
    'linewidth',1.5,'color',cedge1);
if fitWavelength
    hold on;
    ii = linspace(min(xvals),max(xvals),1000);
    lfit = plot(ii,feval(fout,ii));
    legend(lfit,l_text,'interpeter','latex')
end
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
fitPhase=1;
if fitPhase
    
    PhiFit = polyfit(xvals,[stripes.phi]/(2*pi),1);

    phi_fit = PhiFit(1)*xvals+PhiFit(2);
    
    phi_text = sprintf('\\phi =%.5f x + %.2f',PhiFit(1),PhiFit(2));

    end
     
subplot(223);
errorbar(xvals,[stripes.phi]/(2*pi),[stripes.phi_err]/(2*pi),'marker','o',...
    'MarkerFacecolor',cface1,'markeredgecolor',cedge1,'linestyle','none',...
    'linewidth',1.5,'color',cedge1);
if fitPhase
    hold on;
    pfit=plot(xvals,phi_fit);
    legend(pfit,phi_text,'interpeter','latex')
end
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

