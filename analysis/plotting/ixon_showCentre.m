function [hF,hF2]=ixon_showCentre(data,xVar,plt_opts,fit_opts)

if nargin <4 
    fit_opts = struct;
end

if nargin < 3 
    plt_opts = struct;
end

if nargin < 2
    xVar = 'ExecutionDate';
end

if ~isfield(plt_opts,'FigLabel')
    plt_opts.FigLabel = '';
end

doFit = 1;

%% Get Data
Xc = data.Xc;
Yc = data.Yc;

params = [data.Params];
xvals = [params.(xVar)];
PixelSize = data.PixelSize;

px2um = data.PixelSize/data.Magnification;
%% Make Figure

hF=figure('Name',[pad(['iXon ' data.FitType ' centre'],20) plt_opts.FigLabel],...
    'units','pixels','color','w','numbertitle','off');
hF.Position(1)=510;
hF.Position(2)=380;
hF.Position(3)=500;
hF.Position(4)=600;
drawnow;

% Image directory folder string
t=uicontrol('style','text','string',plt_opts.FigLabel,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',6);
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];

uicontrol('style','text','string','PCO','units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',12,'fontweight','bold',...
    'position',[2 2 40 20]);

ixon_resizeFig(hF,t)



%% Track X

if doFit
    hax1=subplot(221);
else
    hax1=subplot(211);
end

set(hax1,'box','on','linewidth',1,'fontsize',10,'units','pixels');
hold on
xlabel([xVar ' (' plt_opts.xUnit ')'],'interpreter','none');
for nn=1:size(Xc,2)
    [cface,cedge] = ixoncolororder(nn);
   plot(xvals,Xc(:,nn),'o','color',cedge,'linewidth',1,'markersize',8,...
       'markerfacecolor',cface,'markeredgecolor',cedge);
end
str='X centre (px)';
text(0.02,.98,str,'units','normalized','fontsize',12,'verticalalignment','cap',...
    'interpreter','latex');

if isequal(xVar,'ExecutionDate')
    datetick('x');
    xlabel('ExecutionDate');
end
%% Table X

if doFit
    a=subplot(222);
    pos=get(a,'position');
    delete(a)

    sTblX=uitable('FontSize',8,'RowName',{},'ColumnName',{},...
        'ColumnEditable',[false false],'units','normalized');
    sTblX.ColumnWidth={100 60};
    sTblX.Position=pos;
    sTblX.Data={[char(0x0394) 'X (px)'],num2str(round(range(Xc(:,nn)),1))};
    drawnow;
end



%% Track Y
if doFit
    hax2=subplot(223);
else
    hax2=subplot(212);
end

set(hax2,'box','on','linewidth',1,'fontsize',10,'units','pixels');
hold on
xlabel([xVar ' (' plt_opts.xUnit ')'],'interpreter','none');
co=get(gca,'colororder');
for nn=1:size(Xc,2)
    [cface,cedge] = ixoncolororder(nn);
   plot(xvals,Yc(:,nn),'o','color',cedge,'linewidth',1,'markersize',8,...
       'markerfacecolor',cface,'markeredgecolor',cedge);
end
str='Y centre (px)';
text(0.02,0.98,str,'units','normalized','fontsize',12,'verticalalignment','cap',...
    'interpreter','latex');

if isequal(xVar,'ExecutionDate')
    datetick('x');
    xlabel('ExecutionDate');
end


%% Table Y
if doFit
    a=subplot(224);
    pos=get(a,'position');
    delete(a)

    sTblY=uitable('FontSize',8,'RowName',{},'ColumnName',{},...
        'ColumnEditable',[false false],'units','normalized');
    sTblY.ColumnWidth={100 60};
    sTblY.Position=pos;
    sTblY.Data={[char(0x0394) 'Y (px)'],num2str(round(range(Yc(:,nn)),1))};
    drawnow;
end

%% Fits
if isfield(ixon_fit_opts,'Center_SineDecay') && ixon_fit_opts.Center_SineDecay && length(xvals)>4
    tVec=linspace(min(xvals),max(xvals),100);    
    
    % X Fit
    axes(hax1);
    fit1=makeSineDecayFit(xvals',Xc(:,nn));
    plot(tVec,feval(fit1,tVec),'r-');  

    sTblX.ColumnWidth={60 60 60};
    sTblY.ColumnWidth={60 60 60};

    cX=coeffvalues(fit1);
    cIntX=confint(fit1);
    
    data={};
    
    %data{1,3}=range(cInt(:,1))/2;
    %data{2,3}=range(cInt(:,2))/2;

    %data{3,3}=1./(range(cInt(:,2))/2);
    %data{4,3}=range(cInt(:,3))/2;
    %data{5,3}=range(cInt(:,4))/2;
    %data{6,3}=range(cInt(:,5))/2;

    
    sTblX.Data={};
    data{1,1}='amp (px)';
    data{2,1}='period';
    data{3,1}='freq';

    data{4,1}='phase (rad)';
    data{5,1}='offset (px)';
    data{6,1}='tau ';

    data{1,2}=cX(1);
    data{2,2}=cX(2);
    data{3,2}=1/cX(2);

    data{4,2}=cX(3);
    data{5,2}=cX(4);
    data{6,2}=cX(5);
    
    data{7,1}='<HTML> &Delta;X (px)</HTML>';
    data{7,2}=range(Xc(:,nn));
    data{8,1}='<HTML> Mean(x) </HTML>';
    data{8,2}=mean(Xc(:,nn));
    
    sTblX.Data=data;
    sTblX.Position(3)=sTblX.Extent(3);
    sTblX.Position(4)=sTblX.Extent(4); 
    
    % Y Fit
    data={};
    
    
    axes(hax2);
    fit2=makeSineDecayFit(xvals',Yc(:,nn));
    cIntY=confint(fit2);
    
    %data{1,3}=range(cInt(:,1))/2;
    %data{2,3}=range(cInt(:,2))/2;

    %data{3,3}=1./(range(cInt(:,2))/2);
    %data{4,3}=range(cInt(:,3))/2;
    %data{5,3}=range(cInt(:,4))/2;
    %data{6,3}=range(cInt(:,5))/2;
    
    
    plot(tVec,feval(fit2,tVec),'r-');  

    cY=coeffvalues(fit2);
    
    sTblY.Data={};
    data{1,1}='amp (px)';
    data{2,1}='period';
    data{3,1}='freq';
    data{4,1}='phase (rad)';
    data{5,1}='offset (px)';
    data{6,1}='tau ';

    data{1,2}=cY(1);
    data{2,2}=cY(2);
    data{3,2}=1/cY(2);

    data{4,2}=cY(3);
    data{5,2}=cY(4);
    data{6,2}=cY(5);
    
    data{7,1}='<HTML> &Delta;Y (px)</HTML>';
    data{7,2}=range(Yc(:,nn));
    data{8,1}='<HTML> Mean(y) </HTML>';
    data{8,2}=mean(Yc(:,nn));
    
    sTblY.Data=data;
    sTblY.Position(3)=sTblY.Extent(3);
    sTblY.Position(4)=sTblY.Extent(4); 
    drawnow;
end

%%

if isfield(ixon_fit_opts,'Center_SineGrow') && ixon_fit_opts.Center_SineGrow && length(xvals)>4
    tVec=linspace(min(xvals),max(xvals),100);    
    
    % X Fit
    axes(hax1);
    fit1=makeSineGrowFit(xvals',Xc(:,nn));
    plot(tVec,feval(fit1,tVec),'r-');  

    sTblX.ColumnWidth={60 60 60};
    sTblY.ColumnWidth={60 60 60};

    cX=coeffvalues(fit1);
    cIntX=confint(fit1);
    
    data={};
    
    %data{1,3}=range(cInt(:,1))/2;
    %data{2,3}=range(cInt(:,2))/2;

    %data{3,3}=1./(range(cInt(:,2))/2);
    %data{4,3}=range(cInt(:,3))/2;
    %data{5,3}=range(cInt(:,4))/2;
    %data{6,3}=range(cInt(:,5))/2;

    
    sTblX.Data={};
    data{1,1}='amp (px)';
    data{2,1}='period';
    data{3,1}='freq';

    data{4,1}='phase (rad)';
    data{5,1}='offset (px)';
    data{6,1}='tau ';

    data{1,2}=cX(1);
    data{2,2}=cX(2);
    data{3,2}=1/cX(2);

    data{4,2}=cX(3);
    data{5,2}=cX(4);
    
    data{7,1}='<HTML> &Delta;X (px)</HTML>';
    data{7,2}=range(Xc(:,nn));
    data{8,1}='<HTML> Mean(x) </HTML>';
    data{8,2}=mean(Xc(:,nn));
    
    sTblX.Data=data;
    sTblX.Position(3)=sTblX.Extent(3);
    sTblX.Position(4)=sTblX.Extent(4); 
    
    % Y Fit
    data={};
    
    
    axes(hax2);
    fit2=makeSineGrowFit(xvals',Yc(:,nn));
    cIntY=confint(fit2);
    
    %data{1,3}=range(cInt(:,1))/2;
    %data{2,3}=range(cInt(:,2))/2;

    %data{3,3}=1./(range(cInt(:,2))/2);
    %data{4,3}=range(cInt(:,3))/2;
    %data{5,3}=range(cInt(:,4))/2;
    %data{6,3}=range(cInt(:,5))/2;
    
    
    plot(tVec,feval(fit2,tVec),'r-');  

    cY=coeffvalues(fit2);
    
    sTblY.Data={};
    data{1,1}='amp (px)';
    data{2,1}='period';
    data{3,1}='freq';
    data{4,1}='phase (rad)';
    data{5,1}='offset (px)';
    data{6,1}='tau ';

    data{1,2}=cY(1);
    data{2,2}=cY(2);
    data{3,2}=1/cY(2);

    data{4,2}=cY(3);
    data{5,2}=cY(4);
    
    data{7,1}='<HTML> &Delta;Y (px)</HTML>';
    data{7,2}=range(Yc(:,nn));
    data{8,1}='<HTML> Mean(y) </HTML>';
    data{8,2}=mean(Yc(:,nn));
    
    sTblY.Data=data;
    sTblY.Position(3)=sTblY.Extent(3);
    sTblY.Position(4)=sTblY.Extent(4); 
    drawnow;
end


%%
if isfield(ixon_fit_opts,'Center_Linear') && ixon_fit_opts.Center_Linear && length(xvals)>1
    tVec=linspace(min(xvals),max(xvals),100);   
    
    D1=Xc(:,nn);    
    D2=Yc(:,nn);

    
    % X Fit
    axes(hax1);
    fit1=polyfit(xvals',D1,1);
    plot(tVec,polyval(fit1,tVec),'r-','linewidth',1);  
    
    sTblX.Data={};
    data{1,1}='slope (px/var)';
    data{2,1}='slope (um/var)';
    data{3,1}='intercept (px)';
    data{4,1}='intercept (um) ';

    pxsize = 16/(80*2);
    
    data{1,2}=fit1(1);
    data{2,2}=fit1(1)*pxsize*1e6;
    data{3,2}=fit1(2);
    data{4,2}=fit1(2)*pxsize*1e6;
    
    data{5,1}='<HTML> &Delta;X (px)</HTML>';
    data{5,2}=range(Xc(:,nn));
    data{6,1}='<HTML> Mean(x) </HTML>';
    data{6,2}=mean(Xc(:,nn));
    
    sTblX.Data=data;
    sTblX.Position(3)=sTblX.Extent(3);
    sTblX.Position(4)=sTblX.Extent(4); 
    
     % Y Fit
    axes(hax2);
    fit2=polyfit(xvals',D2,1);
    plot(tVec,polyval(fit2,tVec),'r-','linewidth',1);  
    
    sTblY.Data={};
    data{1,1}='slope (px/var)';
    data{2,1}='slope (um/var)';
    data{3,1}='intercept (px)';
    data{4,1}='intercept (um) ';

    data{1,2}=fit2(1);
    data{2,2}=fit2(1)*pxsize*1e6;
    data{3,2}=fit2(2);
    data{4,2}=fit2(2)*pxsize*1e6;
    
    data{5,1}='<HTML> &Delta;Y (px)</HTML>';
    data{5,2}=range(Yc(:,nn));
    data{6,1}='<HTML> Mean(y) </HTML>';
    data{6,2}=mean(Yc(:,nn));
    
    sTblY.Data=data;
    sTblY.Position(3)=sTblY.Extent(3);
    sTblY.Position(4)=sTblY.Extent(4); 
    
    

end


% hax1.Position(4)=hax1.Position(4)-15;
% hax2.Position(4)=hax1.Position(4);
%%
ixon_resizeFig(hF2,t)
end

function fitResult=makeSineDecayFit(X,Y,W)

% Guess the amplitude and offset
gA=0.5*range(Y);
gD=(max(Y)+min(Y))*.5;

% Guess the period
iHigh=find((Y-gD)/gA>.8,1);
iLow=find((Y-gD)/gA<-.8,1);
gB=abs(X(iHigh)-X(iLow))*2.2;

gB=8;

minValues=X(Y==min(Y));
maxValues=X(Y==max(Y));
% gB=1*abs(maxValues(1)-minValues(1));
% gB=range(X)/2;




gC=maxValues(1);
gC=pi;
gD=0.5*(max(Y)+min(Y));

gC=pi/2;
gE = range(X);

cosFit=fittype('A*cos(2*pi*t/B+C)*exp(t/E)+D','independent',{'t'},...
    'coefficients',{'A','B','C','D','E'});
options=fitoptions(cosFit);          
        set(options, 'TolFun', 1E-14);
        set(options,'Lower', [0.25*gA,...
            .1*gB,...
            0, ...
            0.75*gD, ...
            0]);
        set(options, 'Upper', [5*gA, ...
            20*gB,...
            2*pi, ...
            1.5*gD, ...
            inf]);            
        set(options, 'StartPoint', [gA, gB,...
            gC,gD, gE]);     
        set(options, 'MaxIter',3000);
        set(options, 'MaxFunEvals',3000);
        set(options,'TolFun',10^-9);
        
        if nargin==3
           set(options,'Weights',W); 
        end        
        
        
        fitResult=fit(X,Y,cosFit,options);      
        
disp(fitResult)
end

function fitResult=makeSineGrowFit(X,Y,W)

% Guess the amplitude and offset
gA=range(Y)/max(X)*0.5;
gD=(max(Y)+min(Y))*.5;

% Guess the period
iHigh=find((Y-gD)/gA>.8,1);
iLow=find((Y-gD)/gA<-.8,1);
gB=abs(X(iHigh)-X(iLow))*2.2;

gB=4;

minValues=X(Y==min(Y));
maxValues=X(Y==max(Y));
% gB=1*abs(maxValues(1)-minValues(1));
% gB=range(X)/2;




gC=maxValues(1);
gC=pi;
gD=0.5*(max(Y)+min(Y));

gC=pi/2;
gE = 0;

cosFit=fittype('(A*t+E)*cos(2*pi*t/B+C)+D','independent',{'t'},...
    'coefficients',{'A','B','C','D','E'});
options=fitoptions(cosFit);          
        set(options, 'TolFun', 1E-14);
        set(options,'Lower', [0.25*gA,...
            .1*gB,...
            0, ...
            0.75*gD,-inf]);
        set(options, 'Upper', [5*gA, ...
            20*gB,...
            2*pi, ...
            1.5*gD,inf]);            
        set(options, 'StartPoint', [gA, gB,...
            gC,gD,gE]);     
        set(options, 'MaxIter',3000);
        set(options, 'MaxFunEvals',3000);
        set(options,'TolFun',10^-9);
        
        if nargin==3
           set(options,'Weights',W); 
        end        
        
        
        fitResult=fit(X,Y,cosFit,options);      
        
disp(fitResult)
end



