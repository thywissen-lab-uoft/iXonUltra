function hF = ixon_showCentre(datain,xVar,plt_opts,fit_opts)

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
Xc = datain.Xc;
Yc = datain.Yc;

params = [datain.Params];
xvals = [params.(xVar)];
PixelSize = datain.PixelSize;

px2um = datain.PixelSize/datain.Magnification;
%% Make Figure

hF=figure('Name',[pad(['iXon ' datain.FitType ' centre'],20) plt_opts.FigLabel],...
    'units','pixels','color','w','numbertitle','off');
hF.Position(1)=510;
hF.Position(2)=380;
hF.Position(3)=700;
hF.Position(4)=500;
drawnow;

% Image directory folder string
t=uicontrol('style','text','string',plt_opts.FigLabel,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',6);
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];

uicontrol('style','text','string',['iXon ' datain.FitType],'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',10,'fontweight','bold',...
    'position',[2 2 80 20]);

ixon_resizeFig(hF,t)



%% Track X

if doFit
    hax1=subplot(221);
else
    hax1=subplot(211);
end

for nn=1:size(Xc,2)
    [cface,cedge] = ixoncolororder(nn);
   plot(xvals,Xc(:,nn),'o','color',cedge,'linewidth',1,'markersize',8,...
       'markerfacecolor',cface,'markeredgecolor',cedge);
   
   
if isequal(xVar,'ExecutionDate')
    datetick('x');
    xlabel('ExecutionDate');
end
   hold on

end
str='X centre (px)';
text(0.02,.98,str,'units','normalized','fontsize',12,'verticalalignment','cap',...
    'interpreter','latex');


set(hax1,'box','on','linewidth',1,'fontsize',10,'units','pixels');
xlabel([xVar ' (' plt_opts.xUnit ')'],'interpreter','none');



%% Table X

if doFit
    a=subplot(222);
    pos=get(a,'position');
    delete(a)

    sTblX=uitable('FontSize',8,'RowName',{},'ColumnName',{},...
        'ColumnEditable',[false false],'units','normalized');
    sTblX.ColumnWidth={100 80};
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


co=get(gca,'colororder');
for nn=1:size(Xc,2)
    [cface,cedge] = ixoncolororder(nn);
   plot(xvals,Yc(:,nn),'o','color',cedge,'linewidth',1,'markersize',8,...
       'markerfacecolor',cface,'markeredgecolor',cedge);
   if isequal(xVar,'ExecutionDate')
    datetick('x');
    xlabel('ExecutionDate');
   end
hold on
end
str='Y centre (px)';
text(0.02,0.98,str,'units','normalized','fontsize',12,'verticalalignment','cap',...
    'interpreter','latex');


set(hax2,'box','on','linewidth',1,'fontsize',10,'units','pixels');
xlabel([xVar ' (' plt_opts.xUnit ')'],'interpreter','none');

%% Table Y
if doFit
    a=subplot(224);
    pos=get(a,'position');
    delete(a)

    sTblY=uitable('FontSize',8,'RowName',{},'ColumnName',{},...
        'ColumnEditable',[false false],'units','normalized');
    sTblY.ColumnWidth={100 80};
    sTblY.Position=pos;
    sTblY.Data={[char(0x0394) 'Y (px)'],num2str(round(range(Yc(:,nn)),1))};
    drawnow;
end

%% Sine 

if isfield(fit_opts,'Center_Sine') && fit_opts.Center_Sine && length(xvals)>4
    tVec=linspace(min(xvals),max(xvals),100);    
    
    % X Fit
    axes(hax1);
    fit1=makeSineFit(xvals',Xc(:,nn));
    plot(tVec,feval(fit1,tVec),'r-');  

    sTblX.ColumnWidth={60 60 60};
    sTblY.ColumnWidth={60 60 60};

    cX=coeffvalues(fit1);
    cIntX=confint(fit1,0.67);
    
    data={};
    
    data{1,3}=range(cIntX(:,1))/2;
    data{2,3}=range(cIntX(:,2))/2;

    data{3,3}=1./(range(cIntX(:,2))/2);
    data{4,3}=range(cIntX(:,3))/2;
    data{5,3}=range(cIntX(:,4))/2;
    data{6,3}=range(cIntX(:,5))/2;

    
    sTblX.Data={};
    data{1,1}='amp (px)';
    data{2,1}='period';
    data{3,1}='freq';

    data{4,1}='phase (rad)';
    data{5,1}='offset (px)';

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
    fit2=makeSineDecayFit(xvals',Yc(:,nn));
    cIntY=confint(fit2,0.67);
    
    data{1,3}=range(cIntY(:,1))/2;
    data{2,3}=range(cIntY(:,2))/2;

    data{3,3}=1./(range(cIntY(:,2))/2);
    data{4,3}=range(cIntY(:,3))/2;
    data{5,3}=range(cIntY(:,4))/2;
    data{6,3}=range(cIntY(:,5))/2;
    
    
    plot(tVec,feval(fit2,tVec),'r-');  

    cY=coeffvalues(fit2);
    
    sTblY.Data={};
    data{1,1}='amp (px)';
    data{2,1}='period';
    data{3,1}='freq';
    data{4,1}='phase (rad)';
    data{5,1}='offset (px)';

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

%% Sine Decay

if isfield(fit_opts,'Center_SineDecay') && fit_opts.Center_SineDecay && length(xvals)>4
    tVec=linspace(min(xvals),max(xvals),100);    
    
    % X Fit
    axes(hax1);
    fit1=makeSineDecayFit(xvals',Xc(:,nn));
    plot(tVec,feval(fit1,tVec),'r-');  

    sTblX.ColumnWidth={60 60 60};
    sTblY.ColumnWidth={60 60 60};

    cX=coeffvalues(fit1);
    cInt=confint(fit1);
    
    data={};
   
    sTblX.Data={};
    data{1,1}='amp (px)';
    data{2,1}='freq';
    data{3,1}='period';
    data{4,1}='phase (rad)';
    data{6,1}='offset (px)';
    data{5,1}='tau ';

    data{1,2}=cX(1);
    data{2,2}=cX(2);
    data{3,2}=1/cX(2);
    data{4,2}=cX(3);
    data{5,2}=cX(4);
    data{6,2}=cX(5);

    data{1,3}=range(cInt(:,1))/2;
    data{2,3}=range(cInt(:,2))/2;
    data{3,3}=(range(cInt(:,2))/2)./(cX(2)).^2;
    data{4,3}=range(cInt(:,3))/2;
    data{5,3}=range(cInt(:,4))/2;
    data{6,3}=range(cInt(:,5))/2;
    
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
    
    plot(tVec,feval(fit2,tVec),'r-');  

    cY=coeffvalues(fit2);
    
    cInt=confint(fit2);
    
    data{1,3}=range(cInt(:,1))/2;
    data{2,3}=range(cInt(:,2))/2;

    data{3,3}=(range(cInt(:,2))/2)./(cY(2)).^2;
    data{4,3}=range(cInt(:,3))/2;
    data{5,3}=range(cInt(:,4))/2;
    data{6,3}=range(cInt(:,5))/2;
    
    sTblY.Data={};
    data{1,1}='amp (px)';
    data{2,1}='freq';
    data{3,1}='period';
    data{4,1}='phase (rad)';
    data{6,1}='offset (px)';
    data{5,1}='tau ';

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

%% Sine Grow

if isfield(fit_opts,'Center_SineGrow') && fit_opts.Center_SineGrow && length(xvals)>4
    tVec=linspace(min(xvals),max(xvals),100);    
    
    % X Fit
    axes(hax1);
    fit1=makeSineGrowFit(xvals',Xc(:,nn));
    plot(tVec,feval(fit1,tVec),'r-');  

    sTblX.ColumnWidth={60 60 60};
    sTblY.ColumnWidth={60 60 60};

    cX=coeffvalues(fit1);
    cIntX=confint(fit1,0.67);
    
    data={};
    
    

    
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
    
    data{1,3}=range(cIntX(:,1))/2;
    data{2,3}=range(cIntX(:,2))/2;

    data{3,3}=(range(cIntX(:,2))/2)/(cX(2)^2);
    data{4,3}=range(cIntX(:,3))/2;
    data{5,3}=range(cIntX(:,4))/2;
    data{6,3}=range(cIntX(:,5))/2;
    
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
    cIntY=confint(fit2,0.67);
    
    
    
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
    
    data{1,3}=range(cIntY(:,1))/2;
    data{2,3}=range(cIntY(:,2))/2;

    data{3,3}=(range(cIntY(:,2))/2)/(cY(2)^2);
    data{4,3}=range(cIntY(:,3))/2;
    data{5,3}=range(cIntY(:,4))/2;
    data{6,3}=range(cIntY(:,5))/2;
    
    
    data{7,1}='<HTML> &Delta;Y (px)</HTML>';
    data{7,2}=range(Yc(:,nn));
    data{8,1}='<HTML> Mean(y) </HTML>';
    data{8,2}=mean(Yc(:,nn));
    
    sTblY.Data=data;
    sTblY.Position(3)=sTblY.Extent(3);
    sTblY.Position(4)=sTblY.Extent(4); 
    drawnow;
end


%% Linear Fit
if isfield(fit_opts,'Center_Linear') && fit_opts.Center_Linear && length(xvals)>1
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

    pxsize = 16/(80);
    
    data{1,2}=fit1(1);
    data{2,2}=fit1(1)*pxsize;
    data{3,2}=fit1(2);
    data{4,2}=fit1(2)*pxsize;
    
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
    data{2,2}=fit2(1)*pxsize;
    data{3,2}=fit2(2);
    data{4,2}=fit2(2)*pxsize;
    
    data{5,1}='<HTML> &Delta;Y (px)</HTML>';
    data{5,2}=range(Yc(:,nn));
    data{6,1}='<HTML> Mean(y) </HTML>';
    data{6,2}=mean(Yc(:,nn));
    
    sTblY.Data=data;
    sTblY.Position(3)=sTblY.Extent(3);
    sTblY.Position(4)=sTblY.Extent(4); 
    
    

end
%%
ixon_resizeFig(hF,t)
end

% function fitResult=makeSineDecayFit(X,Y,W)
% 
% % Guess the amplitude and offset
% gA=0.5*range(Y);
% gD=(max(Y)+min(Y))*.5;
% 
% % Guess the period
% iHigh=find((Y-gD)/gA>.8,1);
% iLow=find((Y-gD)/gA<-.8,1);
% gB=abs(X(iHigh)-X(iLow))*2.2;
% 
% gB=18;
% 
% minValues=X(Y==min(Y));
% maxValues=X(Y==max(Y));
% % gB=1*abs(maxValues(1)-minValues(1));
% % gB=range(X)/2;
% 
% 
% % gA = 20;
% 
% gC=maxValues(1);
% 
% gD=0.5*(max(Y)+min(Y));
% 
% gC=0;
% gC=pi/2;
% gE = range(X)/2;
% % gE = 7;
% 
% cosFit=fittype('A*cos(2*pi*t/B+C)*exp(-t/E)+D','independent',{'t'},...
%     'coefficients',{'A','B','C','D','E'});
% options=fitoptions(cosFit);          
%         % set(options, 'TolFun', 1E-14);
%         % set(options,'Lower', [0.25*gA,...
%         %     .1*gB,...
%         %     0, ...
%         %     0.75*gD, ...
%         %     0]);
%         options.Lower = [0.25*gA,...
%             .1*gB,...
%             0, ...
%             0.75*gD, ...
%             0];
%         % set(options, 'Upper', [5*gA, ...
%         %     20*gB,...
%         %     2*pi, ...
%         %     1.5*gD, ...
%         %     inf]);
%         options.Upper = [5*gA, ...
%             20*gB,...
%             2*pi, ...
%             1.5*gD, ...
%             inf];
%         % set(options, 'StartPoint', [gA, gB,...
%         %     gC,gD, gE]);   
%         options.StartPoint = [gA, gB,...
%             gC,gD, gE];   
%         % set(options, 'MaxIter',3000);
%         options.MaxIter = 3000;
%         % set(options, 'MaxFunEvals',3000);
%         options.MaxFunEvals = 3000;
%         % set(options,'TolFun',10^-9);
%         options.TolFun = 10^-9;
% 
%         if nargin==3
%            set(options,'Weights',W); 
%         end        
% 
%         fitResult=fit(X,Y,cosFit,options);      
% 
% disp(fitResult)
% end

function fitResult=makeSineDecayFit(X,Y,W)

% GUESS : AMPLITUDE
guess_amp = 0.5*range(Y);

% GUESS : OFFSET
guess_off = (max(Y)+min(Y))*.5;

% GUESS : FREQUENCY
dX = diff(sort(unique(X),'ascend'));
dXMin = min(dX);             % Minimum separation sets highest freq
dXMax = max(X) - min(X);     % Total range sets lowest freq

% Set reasonable frequency guess bounds
freq_min = 1.5*(1/dXMax);
freq_max = 0.2*(1/dXMin);
% freq_min = 0.03;
% freq_max = 0.05;
N=1e4;

% Do a correlation measurement with the frequency vector
fVec = linspace(freq_min,freq_max,N);
CC=zeros(N,1);
Yosc = Y(:)-guess_off;
for nn=1:N
    fme = fVec(nn);
    Yf = exp(1i*2*pi*fme*X(:));
    CC(nn) = sum(Yf.*Yosc);
end

% Find the index whose frequency has the highest correlation
[~,ind] = max(abs(CC).^2);

% Get the frequency
guess_freq = fVec(ind);
% Get the phase
guess_phi = atan2(imag(CC(ind)),real(CC(ind))); 

% keyboard

% Tau Guess
guess_tau = max(X) - min(X);
guess_tau = 40;

% In case of sine grow
if max(Y(1:round(length(Y)/2))) < max(Y(round(length(Y)/2):end))
    guess_tau = -guess_tau;
end

% Override Guess
% guess_tau = 1000;           % manual override
% guess_freq = .3; % manual overide


cosFit=fittype('amp*cos(2*pi*freq*t-phi)*exp(-t/tau)+off','independent',{'t'},...
    'coefficients',{'amp','freq','phi','tau','off'});
options=fitoptions(cosFit);          

options.TolFun = 1E-14;
options.Lower  = [...
    0.5*guess_amp,...
    .5*guess_freq,...
    guess_phi-pi, ...
    -abs(guess_tau)*20, ...
    guess_off-100];
options.Upper  = [...
    2*guess_amp, ...
    2.0*guess_freq,...
    guess_phi+pi, ...
    abs(guess_tau)*20, ...
    guess_off+100];
options.StartPoint = [guess_amp, guess_freq,...
    guess_phi,guess_tau, guess_off];
options.MaxIter = 3000;
options.MaxFunEvals = 3000;
options.TolFun = 1E-9;

if nargin==3
   options.Weights = W;
end                

fitResult=fit(X,Y,cosFit,options);             
disp(fitResult);



end




function fitResult=makeSineGrowFit(X,Y,W)

% Guess the amplitude and offset
gA=range(Y)/max(X)*0.5;
gD=(max(Y)+min(Y))*.5;

% Guess the period
iHigh=find((Y-gD)/gA>.8,1);
iLow=find((Y-gD)/gA<-.8,1);
gB=abs(X(iHigh)-X(iLow))*2.2;

gB=20;

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


function fitResult=makeSineFit(X,Y,W)

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

cosFit=fittype('A*cos(2*pi*t/B+C)+D','independent',{'t'},...
    'coefficients',{'A','B','C','D'});
options=fitoptions(cosFit);          
        set(options, 'TolFun', 1E-14);
        set(options,'Lower', [0.25*gA,...
            .1*gB,...
            0, ...
            0.75*gD]);
        set(options, 'Upper', [5*gA, ...
            20*gB,...
            2*pi, ...
            1.5*gD]);            
        set(options, 'StartPoint', [gA, gB,...
            gC,gD]);     
        set(options, 'MaxIter',3000);
        set(options, 'MaxFunEvals',3000);
        set(options,'TolFun',10^-9);
        
        if nargin==3
           set(options,'Weights',W); 
        end        
        
        
        fitResult=fit(X,Y,cosFit,options);      
        
disp(fitResult)
end



