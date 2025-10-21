function hF = ixon_showCounts(data,xVar,plt_opts,fit_opts)

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

%% Sort the data by the parameter given
params=[data.Params];
X=[params.(xVar)]';
xvals = X';
% Natoms = data.N;

% Natoms(:,1) = data.N(:,1) + data.N(:,2);
% Natoms(:,2) = data.N(:,3) + data.N(:,4);
% 
Natoms(:,1) = data.N(:,1);
Natoms(:,2) = data.N(:,2) + data.N(:,3);



opts = fit_opts;

%% Fitting

%% Exponential Decay Fit

if isfield(opts,'NumberExpFit') && opts.NumberExpFit && size(Natoms,1)>2
    myfit=fittype('A*exp(-t/tau)','coefficients',{'A','tau'},...
    'independent','t');
    opt=fitoptions(myfit);
    
    % Get some initial guesses
    tau0=max(X)/2;
%     tau0=0.1;
    
    fout_exp={};
    for nn=1:size(Natoms,2)  
        A0=max(Natoms(:,nn));
        
        % Assign start point
        opt.StartPoint=[A0 tau0];
        fout_exp{nn}=fit(X,Natoms(:,nn),myfit,opt);
    end
end
% 
% if isfield(opts,'NumberExpOffsetFit') && opts.NumberExpOffsetFit && size(Natoms,1)>3
%     myfit=fittype('A*exp(-t/tau)+B','coefficients',{'A','tau','B'},...
%     'independent','t');
%     opt=fitoptions(myfit);
%     
%     fout_exp={};
%     
%     for nn=1:size(Natoms,2)
%        A0=range(Natoms(:,nn));
%        B0=min(Natoms(:,nn));
%        tau0=range(X)/4;
%        
%        opt.StartPoint=[A0 tau0 B0];
%        opt.Lower=[0 0 0];
%        fout_exp{nn}=fit(X,Natoms(:,nn),myfit,opt);       
%     end       
% end

%% Lorentzian Fits
% Addding lorentzian fits, find a better way to put this in the code in
% the future, future me. Thanks <3

% doLorentzianFit=opts.NumberLorentzianFit;
% fouts_lorentz={};
% if doLorentzianFit
%     for rr=1:size(Natoms,2)
%         myfit=fittype('A*((x-x0).^2+(G/2).^2).^(-1).*(G/2)^2','coefficients',{'A','G','x0'},...
%             'independent','x');
%         opt=fitoptions(myfit);
%         A0=max(Natoms(:,rr));
%         G0=range(X)/2;
% 
%         inds=[Natoms(:,rr)>.8*max(Natoms(:,rr))];
%         x0=mean(X(inds));       
%         opt.StartPoint=[A0 G0 x0];    
%         fout_lorentz=fit(X,Natoms(:,rr),myfit,opt);
%         fouts_lorentz{rr}=fout_lorentz;
%     end
% end
% 


%% Make Figure
hF=figure('Name',[pad(['ixon ' data.FitType ' number'],20) plt_opts.FigLabel],...
    'units','pixels','color','w','Menubar','figure','Resize','on',...
    'numbertitle','off');
hF.Position=[5 380 500 300];clf;

% Image directory folder string
t=uicontrol('style','text','string',plt_opts.FigLabel,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',6);
drawnow;
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];

% Draw iXon label
uicontrol('style','text','string',['iXon ' data.FitType],'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',10,'fontweight','bold',...
    'position',[2 2 80 20]);

% Make axis
hax=axes;
set(hax,'box','on','linewidth',1,'fontsize',10,'xgrid','on','ygrid','on');
hold on
xlabel([xVar ' (' plt_opts.xUnit ')'],'interpreter','none');
ylabel([data.FitType ' number']);
        
% Plot the data
for nn=1:size(Natoms,2)

    % if nn==1
    %         yyaxis left;
    %         ylabel([data.FitType ' number 1']);
    % else
    %         yyaxis right;
    %         ylabel([data.FitType ' number p']);
    % end

    [cface,cedge] = ixoncolororder(nn);
   plot(X,Natoms(:,nn),'o','color',cedge,'linewidth',1,'markersize',8,...
       'markerfacecolor',cface,'markeredgecolor',cedge);
end
% 
if (opts.NumberExpFit || opts.NumberExpOffsetFit) && size(Natoms,1)>3
    strs={};
    xx=linspace(0,max(X),1000);
    
    for nn=1:size(Natoms,2)

        
        [cface,cedge] = ixoncolororder(nn);
        
        pExp(nn)=plot(xx,feval(fout_exp{nn},xx),'-','linewidth',1,...
            'color',0.8*cedge); 
        
        if opts.NumberExpFit        
            fstr=['$N_0 = ' num2str(round(fout_exp{nn}.A),'%.2e') '$' newline ...
                '$\tau = ' num2str(round(fout_exp{nn}.tau,2),'%.2e') ' $'];
        else
            fstr=['$\Delta N = ' num2str(round(fout_exp{nn}.A),'%.2e') '$' newline ...
                '$\tau = ' num2str(round(fout_exp{nn}.tau,1),'%.2e') ' $' newline ...
                '$N_0 = ' num2str(round(fout_exp{nn}.B),'%.2e') '$'];
        end
        strs{nn}=fstr;
    end       
    legend(pExp,strs,'interpreter','latex','location','best','fontsize',8);
    hax.YLim(1)=0;
end
% 
% if doLorentzianFit
%     xx=linspace(min(X),max(X),100);
%     legStr={};
%     
%     for rr=1:length(fouts_lorentz)
%         fout_lorentz=fouts_lorentz{rr};
%         pFs(rr)=plot(xx,feval(fout_lorentz,xx),'-','linewidth',2,'color',...
%         co(rr,:)*.8);
% 
%         str=['$N_0 = ' num2str(round(fout_lorentz.A),'%.2e') '$' newline ...
%             '$\mathrm{FWHM} = ' num2str(round(fout_lorentz.G,3)) ' $' newline ...
%             '$x_0 = ' num2str(round(fout_lorentz.x0,3)) '$'];
%         legStr{rr}=str;
%     end
%     legend(pFs,legStr,'interpreter','latex','location','best','fontsize',8);
% 
%     hax.YLim(1)=0;
% end

if isequal(xVar,'ExecutionDate')
    datetick('x');
    xlabel('ExecutionDate');
end
yl=get(gca,'YLim');
set(gca,'YLim',[0 yl(2)]);


if size(Natoms,2)==2    
    yyaxis right
    set(gca,'YColor',[.5 .5 .5]);
    
    
    [XU, ~, subs] = unique(X);
    
    R = accumarray(subs, Natoms(:,2)./Natoms(:,1), [], @mean);
    Rerr = accumarray(subs, Natoms(:,2)./Natoms(:,1), [], @std);


       errorbar(XU,R,Rerr,'o','color','k','linewidth',1,'markersize',8,...
       'markerfacecolor',[.7 .7 .7],'markeredgecolor','k');
   ylabel('ratio');
   yyaxis right
   
end
ixon_resizeFig(hF,t,[hax]);

%% Sine Decay

if isfield(fit_opts,'Number_SineDecay') && fit_opts.Number_SineDecay && length(xvals)>4

    yyaxis left

    tVec=linspace(min(xvals),max(xvals),100);   
    
    % X Fit
    axes(hax);
    fit1=makeSineDecayFit(xvals',Natoms(:,nn));
    plot(tVec,feval(fit1,tVec),'r-');  

    % sTblX.ColumnWidth={60 60 60};
    % sTblY.ColumnWidth={60 60 60};
    % 
    % cX=coeffvalues(fit1);
    % cInt=confint(fit1);
    % 
    % data={};
    % 
    % sTblX.Data={};
    % data{1,1}='amp (px)';
    % data{2,1}='freq';
    % data{3,1}='period';
    % data{4,1}='phase (rad)';
    % data{6,1}='offset (px)';
    % data{5,1}='tau ';
    % 
    % data{1,2}=cX(1);
    % data{2,2}=cX(2);
    % data{3,2}=1/cX(2);
    % data{4,2}=cX(3);
    % data{5,2}=cX(4);
    % data{6,2}=cX(5);
    % 
    % data{1,3}=range(cInt(:,1))/2;
    % data{2,3}=range(cInt(:,2))/2;
    % data{3,3}=(range(cInt(:,2))/2)./(cX(2)).^2;
    % data{4,3}=range(cInt(:,3))/2;
    % data{5,3}=range(cInt(:,4))/2;
    % data{6,3}=range(cInt(:,5))/2;
    % 
    % data{7,1}='<HTML> &Delta;X (px)</HTML>';
    % data{7,2}=range(Xc(:,nn));
    % data{8,1}='<HTML> Mean(x) </HTML>';
    % data{8,2}=mean(Xc(:,nn));
    % 
    % sTblX.Data=data;
    % sTblX.Position(3)=sTblX.Extent(3);
    % sTblX.Position(4)=sTblX.Extent(4); 
    % 
    % % Y Fit
    % data={};
    % 
    % 
    % axes(hax2);
    % fit2=makeSineDecayFit(xvals',Yc(:,nn));
    % 
    % plot(tVec,feval(fit2,tVec),'r-');  
    % 
    % cY=coeffvalues(fit2);
    % 
    % cInt=confint(fit2);
    % 
    % data{1,3}=range(cInt(:,1))/2;
    % data{2,3}=range(cInt(:,2))/2;
    % 
    % data{3,3}=(range(cInt(:,2))/2)./(cY(2)).^2;
    % data{4,3}=range(cInt(:,3))/2;
    % data{5,3}=range(cInt(:,4))/2;
    % data{6,3}=range(cInt(:,5))/2;
    % 
    % sTblY.Data={};
    % data{1,1}='amp (px)';
    % data{2,1}='freq';
    % data{3,1}='period';
    % data{4,1}='phase (rad)';
    % data{6,1}='offset (px)';
    % data{5,1}='tau ';
    % 
    % data{1,2}=cY(1);
    % data{2,2}=cY(2);
    % data{3,2}=1/cY(2);
    % 
    % data{4,2}=cY(3);
    % data{5,2}=cY(4);
    % data{6,2}=cY(5);
    % 
    % data{7,1}='<HTML> &Delta;Y (px)</HTML>';
    % data{7,2}=range(Yc(:,nn));
    % data{8,1}='<HTML> Mean(y) </HTML>';
    % data{8,2}=mean(Yc(:,nn));
    % 
    % sTblY.Data=data;
    % sTblY.Position(3)=sTblY.Extent(3);
    % sTblY.Position(4)=sTblY.Extent(4); 
    % drawnow;
end



end

%% Fit funcs

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

% guess_freq = 0.17;

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
