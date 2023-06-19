%% pco_analysis_custom.m
%
% This script runs customized analysis on the box, gaussian, or erf fits.
% In particular, it enables you to customized how the analysis is
% performed.

%% Custom Analysis Data Source
% Select the data source

data = ixondata

%% Fit Flags

%%%%%%%%%%%%%%%% Fit Flags
T2exp=0;
Rabi_oscillation = 0;

gauss_single=0;
gauss_4=0;
gauss_neg_double=0;
gauss_neg_single=0;
gauss_double = 0;
expdecay = 1;

lorentz_neg_single=0;    
lorentz_neg_double=0;  

lorentz_single=0;
lorentz_double=0;    
lorentz_triple=0;    

lorentz_asym_single= 0;
lorentz_asym_double= 0;

fit_lorentz_assymetric_4=0;

UniPowerLaw = 0;

%% Generate Custom Data

doCustomX = 1;

process_data = struct;
process_data.Source = data;
process_data.FitType = 'custom';
process_data.Ratio_79=1;    

% Scale the atom number in each box if necessary
N = data.BoxCount;
process_data.Natoms = N;     

% Get the default X data;
process_data.X = [data.Params.(ixon_xVar)];
process_data.XLabel = ixon_xVar;
process_data.XUnit = data.Units(1);

if doCustomX
    if isfield(data.Params(1),'HF_FeshValue_Initial_lattice') && ...
            isfield(data.Params(1),'HF_zshim_Initial_Lattice')
        Bfb   = data.Params(1).HF_FeshValue_Initial_Lattice;
        Bshim = data.Params(1).HF_zshim_Initial_Lattice*2.35;
        Boff  = 0.11;

        B = Bfb + Bshim + Boff;

        % Choose the mf States
        mF1 = -7/2;
        mF2 = -9/2;
        x0 = abs((BreitRabiK(B,9/2,mF1)-BreitRabiK(B,9/2,mF2)))/6.6260755e-34/1E6; 
    end
    
        B = 197 + 0 + 0.11;

        % Choose the mf States
        mF1 = -7/2;
        mF2 = -9/2;
        x0 = abs((BreitRabiK(B,9/2,mF1)-BreitRabiK(B,9/2,mF2)))/6.6260755e-34/1E6; 
        disp(x0)

    switch pco_xVar
        case 'Raman_AOM3_freq'
            X=process_data.X;
            X = 2*X - 80;
            X = X - x0;   
            X = X*1e3;
            xstr=['frequency - ' num2str(round(abs(x0),4))  ' MHz (kHz)']; 
            xunit = 'kHz';
        case 'pulse_time' %'Pulse_Time'
            X=process_data.X;
            xstr='pulse time (ms)';    
            xunit = 'ms';
       case 'rf_rabi_time_HF'
            X=process_data.X;
            xstr='pulse time (ms)';    
            xunit = 'ms';
        case 'rf_freq_HF'
            X=process_data.X;
            X = X - x0;   
            X = X*1e3;
            xstr=['frequency - ' num2str(round(abs(x0),4))  ' MHz (kHz)']; 
            xunit = 'kHz';
        case 'rf_tof_freq'
          X=process_data.X;
            X = X - x0;   
            X = X*1e3;
            xstr=['frequency - ' num2str(round(abs(x0),4))  ' MHz (kHz)'];  
            xunit = 'kHz';
        otherwise
            X = process_data.X;
            xstr = pco_xVar;
            xunit = '';
    end 
    process_data.X = X;
    process_data.XLabel = xstr;        
    process_data.XUnit = xunit;
end    


%% Rabi Oscillations

rabi_opts=struct;
rabi_opts.FigLabel  = FigLabel;
rabi_opts.xUnit     = pco_unit;

rabi_opts.Ratio_79=0.62;
rabi_opts.Ratio_79=1;

rabi_opts.Guess=[.9 13 1]; % [probability transfer, freq, t2 time]
rabi_opts.Sign='auto'; % Automatic fit sign

T = process_data.X;
N = process_data.Natoms;

if doRabiContrast && size(data.Natoms,1)>4
    
    % Calculate the contrast for two boxes
    if size(data.Natoms,2)==2   
        C = (N(:,1)-N(:,2))./(N(:,1)+N(:,2));
        Ntot = N(:,1)+N(:,2);
    end
    
    % Calculate the contrast for three boxes
    if size(data.Natoms,3)==3
        C = (N(:,1)+N(:,2)-N(:,3))./(N(:,1)+N(:,2)+N(:,3));
        Ntot = N(:,1)+N(:,2)+N(:,3);
    end  
    
    % Add additional lines for different counting schemes
    [hF_rabi_contrast,rabi_contrast]=rabiOscillationsContrast(data,pco_xVar,rabi_opts);
    if doSave;saveFigure(hF_rabi_contrast,[data_source '_rabi_oscillate_contrast'],saveOpts);end
end

if doRabiAbsolute && size(data.Natoms,1)>4    
    [hF_rabi_raw,rabi_absolute]=rabiOscillationsAbsolute(data,pco_xVar,rabi_opts);
    if doSave;saveFigure(hF_rabi_raw,[data_source '_rabi_oscillate_raw'],saveOpts);end    
end
    
%% Process X Data

% Center frequency for expected RF field (if relevant)
% Calibrated 2021/09/25-26
% Bfb   = data.Params(1).HF_FeshValue_Initial_Lattice;
% Bfb   = data.Params(1).HF_FeshValue_Initial_ODT;
Bshim =0;
% Bfb   = data.Params(1).HF_FeshValue_Spectroscopy;
% Bshim = data.Params(1).HF_zshim_Initial_Lattice*2.35;

Boff  = 0.11;

% B = Bfb + Bshim + Boff;
 
B= 20 +0.11;

% Choose the mf States
mF1 = -7/2;
mF2 = -9/2;

x0 = abs((BreitRabiK(B,9/2,mF1)-BreitRabiK(B,9/2,mF2)))/6.6260755e-34/1E6; 

switch pco_xVar
    case 'Raman_AOM3_freq'
        % Calculate relative to the Raman condition
        X=data.X;
        X = 2*X - 80;  %Raman AOM condition
        X = X - x0;   
        X = X*1e3;
        xstr=['frequency - ' num2str(round(abs(x0),4))  ' MHz (kHz)']; 
    case 'Pulse_Time'
        X=data.X;
        xstr=['pulse time (ms)'];    
   case 'rf_rabi_time_HF'
        X=data.X;
        xstr=['pulse time (ms)'];    
    case 'rf_freq_HF'
        X=data.X;
        X = X - x0;   
        X = X*1e3;
        xstr=['frequency - ' num2str(round(abs(x0),4))  ' MHz (kHz)'];        
    case 'rf_rabi_freq_HF'
        X=data.X;
        X = X - x0;   
        X = X*1e3;
        xstr=['frequency - ' num2str(round(abs(x0),4))  ' MHz (kHz)']; 
    case 'rf_tof_freq'
      X=data.X;
        X = X - x0;   
        X = X*1e3;
        xstr=['frequency - ' num2str(round(abs(x0),4))  ' MHz (kHz)'];  
    otherwise
        X = data.X;
        xstr = pco_xVar;
end    

%% Define the Y Data
    
    % Scale atom number if it imaged at -7
    Ratio_79=1;0.9;
    N = data.Boxcount;
    for nn=1:size(data.Natoms,2)
       if mean(data.Yc(:,nn))>1024
           N(:,nn) = N(:,nn)/Ratio_79;
       end        
    end
    
    % Default total atom number is just the sum
    Ntot = sum(N,2);     
     N(:,2) = N(:,2)*8.6/7.2; %fudge factor

%      dataMode= 2;      
%      
%      Ytype = 'relative';
% 
%      switch dataMode
%          case 0     
%              Y=(N(:,1)-N(:,2))./N(:,1);
%              ystr=['\Delta N97/N9'];
%              figName='custom';
%          case 1
%             Y = N(:,2);
%             ystr=['N_7'];
%             figName=ystr;
%             Ytype='absolute';
%          case 2     
%             Y= N(:,1);
%             ystr=['N_9'];
%             figName=ystr;
%             Ytype='absolute';
%          case 3
%             Y=N(:,1)+N(:,2);
%             ystr=['N_9+N_7'];
%             figName=ystr;
%             Ytype='absolute';
%          case 4
%             Y=N(:,2)./(N(:,1)+N(:,2));
% %             ystr=['Transfer Fraction'];
%             ystr=['N_9/(N_7+N_9)'];
%             figName='Transfer Fraction';
%          case 5 % random customized stuffs 
%              Y =(N(:,1)- N(:,2))./(N(:,1)+N(:,2));
%              ystr=['(N_9-N_7)/(N_7+N_9)'];
%              figName='custom';
%          case 6
%              Y=N(:,2)./N(:,1);
%              ystr=['N_7/N_9'];
%              figName='79 ratio';             
%           case 7
%              Y=N(:,1)./N(:,2);
%              ystr=['N_9/N_7'];
%              figName='79 ratio';             
%          case 8
%              Y =(N(:,1)+N(:,2))./(N(:,1)+N(:,2)+N(:,3));
%              ystr=['y excited fraction'];
%              figName='y excited ratio';             
%           case 9 % 
%              Y =(N(:,2)- N(:,1))./(N(:,1)+N(:,2));
%              ystr=['(N_7-N_9)/(N_7+N_9)'];
%              figName='custom7';             
%           case 10 
%              Y =(N(:,2))./(N(:,1)+N(:,2));
%              ystr=['(N_7)/(N_7+N_9)'];
%              figName='custom72';
%      end

    Y = data.BoxCount;

    custom_data=struct;    
    custom_data.Source = data;
    custom_data.X=X;
    custom_data.Xstr=xstr;   
    custom_data.Ratio_79=Ratio_79;     
    custom_data.N = N;
    custom_data.Y = Y;
    custom_data.YStr = ystr;
    
%% Plot the raw Data

% Get the unique values and error bars
[ux,ia,ib]=unique(X);    
Yu=zeros(length(ux),2);    
for kk=1:length(ux)
    inds=find(X==ux(kk));
    Yu(kk,1)=mean(Y(inds));
    Yu(kk,2)=std(Y(inds));       
end

hFB=figure;
hFB.Color='w';
hFB.Name='box custom';

hFB.Name=figName;
hFB.Position=[400 400 600 350];

co=get(gca,'colororder');    


% Image directory folder string
t=uicontrol('style','text','string',FigLabel,'units','pixels',...
    'backgroundcolor','w','horizontalalignment','left','fontsize',6);
t.Position(3:4)=[hFB.Position(3) 7];
t.Position(1:2)=[5 hFB.Position(4)-t.Position(4)];




% Plot the data with error bars
errorbar(ux,Yu(:,1),Yu(:,2),'o','markerfacecolor',co(1,:),...
    'markeredgecolor',co(1,:)*.5,...
    'linewidth',2,'markersize',8);    

xlabel(xstr,'interpreter','latex');    
ylabel(ystr);

set(gca,'fontsize',10,'linewidth',1,'box','on',...
    'xgrid','on','ygrid','on');
yL=get(gca,'YLim');

hold on    
xlim([min(X) max(X)]);

if isequal(Ytype,'absolute')
    ylim([0 yL(2)]);
end

%% Exponential Fit
    if T2exp
        myfit=fittype('A+(1-A)*exp(-pi*t/tau)',...
            'coefficients',{'A','tau'},...
            'independent','t');
        
        myfit=fittype('0.5+0.5*exp(-pi*t/tau)-A',...
            'coefficients',{'A','tau'},...
            'independent','t');
         
        % Fit options and guess
        opt=fitoptions(myfit);        
        Ag = 0.5;
        taug = median(X);
        G=[Ag taug];        
        opt.StartPoint=G;
  
        % Perform the fit
        fout=fit(X,Y,myfit,opt)
        
        % Plot the fit
        tt=linspace(0,max(X),1000);
        xlim([0 max(X)]);
        pF=plot(tt,feval(fout,tt),'r-','linewidth',1);
        lStr=['$ \tau = ' num2str(round(fout.tau,3)) '~\mathrm{ms}$'];
        legend(pF,lStr,'location','best','interpreter','latex');        
        str = '$A+(1-A)\exp(-\pi t/\tau)$';
        t=text(.02,.03,str,'units','normalized',...
            'fontsize',10,'interpreter','latex');
    end
   
    %% Negative Double Gauss
    if length(X)>8 && gauss_neg_double
        myfit=fittype(['bg-A1*exp(-(x-x1).^2/(2*s1.^2))- ' ...
            'A2*exp(-(x-x2).^2/(2*s2^2))'],...
            'coefficients',{'A1','s1','x1','A2','s2','x2','bg'},...
            'independent','x');
        opt=fitoptions(myfit);
        % Background is max
        bg=max(Y);
        % Find center
        [Ymin,ind]=min(Y);
        A=bg-Ymin;
        xC=X(ind);
        
        % Assign guess        
        xC1 = 0;
        xC2 = 40;
        G=[A 15 xC1 A/2 50 xC2 bg];
        
        opt.StartPoint=G;
        opt.Robust='bisquare';
%         opt.Lower=[0 0 -inf 0 0 -inf 0];
        % Perform the fit
        fout=fit(X,Y,myfit,opt);
        disp(fout);
        ci = confint(fout,0.95);
        disp(ci)
        % Plot the fit
        tt=linspace(min(X),max(X),1000);
        pF=plot(tt,feval(fout,tt),'r-','linewidth',1);
%         ylim([-0.1 2])
        lStr=['xC=(' num2str(round(fout.x1,2)) ' � ' num2str(abs(round(ci(1,3)-fout.x1,2))) ', '...
            num2str(round(fout.x2,2)) ' � ' num2str(abs(round(ci(1,6)-fout.x2,2))) ')' ...
            ' \sigma=(' num2str(round(fout.s1,1)) ', ' num2str(round(fout.s2,1)) ')' ];
        legend(pF,lStr,'location','best');
    end
    
 %% Negative Gauss
 
    if length(X)>4 && gauss_neg_single
        myfit=fittype('bg-A1*exp(-(x-x1).^2/G1.^2)',...
            'coefficients',{'A1','G1','x1','bg'},'independent','x');
        opt=fitoptions(myfit);
        % Background is max
        bg=max(Y);
        % Find center
        [Ymin,ind]=min(Y);
        A=bg-Ymin;
        xC=X(ind);
        % Assign guess
        G=[A 10 38 bg];
        opt.StartPoint=G;
        opt.Robust='bisquare';
        opt.Lower=[0 0 -inf 0 0 -inf 0];
        % Perform the fit
        fout=fit(X,Y,myfit,opt);
        disp(fout);
        ci = confint(fout,0.95);
        disp(ci)
        % Plot the fit
        tt=linspace(min(X),max(X),1000);
        pF=plot(tt,feval(fout,tt),'r-','linewidth',1);
        lStr=['xC=(' num2str(round(fout.x1,2)) '�' num2str(abs(round(ci(1,3)-fout.x1,2))) ')'',' ...
            ' FWHM=(' num2str(round(fout.G1,1)) ')'','...
            ' A=(' num2str(round(fout.A1,1)) ')'];
        legend(pF,lStr,'location','best');
    end
    
%% Negative Lorentz Double
    
    if length(X)>4 && lorentz_neg_double
        X= reshape(X,length(X),1);
        myfit=fittype('bg-A1*(G1/2).^2*((x-x1).^2+(G1/2).^2).^(-1)-A2*(G2/2).^2*((x-x2).^2+(G2/2).^2).^(-1)',...
            'coefficients',{'A1','G1','x1','A2','G2','x2','bg'},...
            'independent','x');
        opt=fitoptions(myfit);
        
        % Background is max
        bg=max(Y);
        
        % Find center
        [Ymin,ind]=min(Y);
        A=bg-Ymin;        
        xC=X(ind);
        
        % Assign guess
        xC1 = -43.5;
        xC2 = -43.7;
%         G=[A/2 20 xC1 A 20 xC2 bg];        
        G=[A/2 0.2 xC1 A 0.2 xC2 bg];
        
        opt.StartPoint=G;
        opt.Robust='bisquare';
        opt.Lower=[0 0 -44.18 0 0 -inf 0];
        
        % Perform the fit
        fout=fit(X,Y,myfit,opt);
        disp(fout);

        % Plot the fit
        tt=linspace(min(X),max(X),1000);
        pF=plot(tt,feval(fout,tt),'r-','linewidth',1);
        lStr=['xC=(' num2str(round(fout.x1,2)) ',' num2str(round(fout.x2,2)) ')' ...
            ' FWHM=(' num2str(round(fout.G1,2)) ',' num2str(round(fout.G2,2)) ')' ];
        legend(pF,lStr,'location','best');
        
        custom_data.Fit=fout;
    end
    
%% Gauss
    if length(X)>10 && gauss_4
        y=@(x,A,s,x0) A*exp(-(x-x0).^2./(2*s^2));
        
        yTot = @(a1,a2,a3,a4,...
                s1,s2,s3,s4,...
                x1,x2,x3,x4,...
                bg,x) ...
            y(x,a1,s1,x1) + ...
            y(x,a2,s2,x2) + ...
            y(x,a3,s3,x3) + ...
            y(x,a4,s4,x4) + bg;           
        
        myfit=fittype(@(a1,a2,a3,a4,...
                s1,s2,s3,s4,...
                x1,x2,x3,x4,...
                bg,x) yTot(a1,a2,a3,a4,...
                s1,s2,s3,s4,...
                x1,x2,x3,x4,...
                bg,x),'coefficients',{'a1','a2','a3','a4',...
                's1','s2','s3','s4',...
                'x1','x2','x3','x4',...
                'bg'},...
            'independent','x'); 
        
        opt=fitoptions(myfit);
        
        opt.StartPoint = zeros(1,13);
        
        % ampitude
        opt.StartPoint(1) = -2.5e4;
        opt.StartPoint(2) = -1.2e4;
        opt.StartPoint(3) = -1e4;
        opt.StartPoint(4) = -2e4;
        
        % sigma
        opt.StartPoint(5) = 10;
        opt.StartPoint(6) = 10;
        opt.StartPoint(7) = 10;
        opt.StartPoint(8) = 10;
        
        % center
        opt.StartPoint(9) = -4.5;
        opt.StartPoint(10) = 60;
        opt.StartPoint(11) = 126;
        opt.StartPoint(12) = 190;
        
        % bkacground
        opt.StartPoint(end) = 4E4;        
        
        fout=fit(X,Y,myfit,opt);
%         ci = confint(fout_lorentz,0.95);   
        
        XF=linspace(min(X)-5,max(X)+5,1000);
        xlim([min(X)-0.1 max(X)+0.1]);
        pExp=plot(XF,feval(fout,XF),'r-','linewidth',2);
%         str=['$f_0 = ' num2str(round(fout_lorentz.x0,2)) '\pm' num2str(round((ci(2,2)-ci(1,2))/2,2)) '$ kHz' newline ...
%             '$\mathrm{FWHM} = ' num2str(round(abs(fout_lorentz.G),2)) ' $ kHz'];
%         legend(pExp,{str},'interpreter','latex','location','best','fontsize',8);         
    end
    
 %% Double Gauss
 
    if length(X)>4 && gauss_double
        myfit=fittype('bg+A1*exp(-(x-x1).^2/G1.^2)+A2*exp(-(x-x2).^2/G2.^2)',...
            'coefficients',{'A1','G1','x1','A2','G2','x2','bg'},'independent','x');
        opt=fitoptions(myfit);
        % Background is max
        bg=min(Y);
        % Find center
        [Ymin,ind]=min(Y);
        A=bg-Ymin;
        xC=X(ind);
        % Assign guess
        G=[A 10 -3 A/5 10 30 bg];
        opt.StartPoint=G;
        opt.Robust='bisquare';
%         opt.Lower=[0 0 -inf 0 0 -inf 0];
        % Perform the fit
        fout=fit(X,Y,myfit,opt);
        disp(fout);
        ci = confint(fout,0.95);
        disp(ci)
        % Plot the fit
        tt=linspace(min(X),max(X),1000);
        pF=plot(tt,feval(fout,tt),'r-','linewidth',1);
        
        str=['$(A_i,f_i,\Gamma_i)$' newline ....
            '$(' num2str(fout.A1,2) ',' ...
            num2str(round(fout.x1,1)) '\pm ' ...
            num2str(round((ci(2,3)-ci(1,3))/2,1)) ...
            ',' num2str(round(fout.G1,1)) ')$' newline ...
            '$(' num2str(fout.A2,2) ',' ...
            num2str(round(fout.x2,1)) '\pm ' ...
            num2str(round((ci(2,6)-ci(1,6))/2,1)) ...
            ',' num2str(round(fout.G2,1)) ')$'];
            
        
%         str=['$f_1 = ' num2str(round(fout.x1,2)) '\pm ' num2str(round((ci(2,3)-ci(1,3))/2,2)) '$ kHz' newline ...
%             '$\mathrm{FWHM} = ' num2str(round(abs(fout.G1),2)) ' $ kHz' newline ...
%             '$f_2 = ' num2str(round(fout.x2,2)) '\pm ' num2str(round((ci(2,6)-ci(1,6))/2,2)) '$ kHz' newline ...
%             '$\mathrm{FWHM} = ' num2str(round(abs(fout.G2),2)) ' $ kHz' newline...
%             'A1 =' num2str(round(fout.A1,2)) newline...
%             'A2 =' num2str(round(fout.A2,2))];
        legend(pF,str,'location','best','interpreter','latex');
        
    end
%% Exponential decay

if length(X)>4 && expdecay
    
    myfit=fittype('A0 + A1*exp(-1*t/tau)',...
        'coefficients',{'A0','A1','tau'},...
        'independent','t');

    % Fit options and guess
    opt=fitoptions(myfit);        
    Ag = max(Y);
    A0 = min(Y);
    taug = median(X);
    G=[A0 Ag taug];        
    opt.StartPoint=G;
    opt.Lower = [0 0 0];
    opt.Upper = [0 inf inf];

    % Perform the fit
    fout=fit(X,Y,myfit,opt);

    % Plot the fit
    tt=linspace(0,max(X),1000);
    xlim([0 max(X)]);
    pF=plot(tt,feval(fout,tt),'r-','linewidth',1);
    lStr=['$ \tau = ' num2str(round(fout.tau,3)) '~\mathrm{ms}$' ...
        newline ...
        '$A_0 = ' num2str(round(fout.A0,3),'%.3e') '$' ... 
        newline ...
        '$A_1 = ' num2str(round(fout.A1,3),'%.3e') '$'];
    legend(pF,lStr,'location','best','interpreter','latex');        
    str = '$A_0+ A_1\exp(-t/\tau)$';
    t=text(.02,.03,str,'units','normalized',...
        'fontsize',10,'interpreter','latex');
end

    %% Negative Lorentzian
    if length(X)>4 && lorentz_neg_single
        myfit=fittype('bg-A*(G/2).^2*((x-x0).^2+(G/2).^2).^(-1)',...
            'coefficients',{'A','G','x0','bg'},...
            'independent','x');
        opt=fitoptions(myfit);
        
        % Background is max
        bg=max(Y);
        
        % Find center
        [Ymin,ind]=min(Y);
        A=bg-Ymin;   
        A=range(Y);
        
%         A=4e4;
%         bg= 5e4;
        xC=X(ind);
%         xlim([-49.7 -49.45]);
%         xC=-41.6;

        % Assign guess
        G=[A 0.025 xC bg];
        opt.StartPoint=G;

        % Perform the fit
        fout=fit(X,Y,myfit,opt);
        disp(fout);

        % Plot the fit
        tt=linspace(min(X),max(X),1000);
        pF=plot(tt,feval(fout,tt),'r-','linewidth',1);
          lStr=['xC=(' num2str(round(fout.x0,4),'%.4f') ')' ...
            ' FWHM=(' num2str(round(fout.G,4),'%.4f') ')' ];
        legend(pF,lStr,'location','best');
        
        cc= confint(fout);
        dG = abs((cc(2,2)-cc(1,2))/2);
        disp(dG)
    end
    
    %% Double Loretnzian
    if length(X)>4 && lorentz_double
        myfit=fittype('bg+A1*(G1/2).^2*((x-x1).^2+(G1/2).^2).^(-1)+A2*(G2/2).^2*((x-x2).^2+(G2/2).^2).^(-1)',...
            'coefficients',{'A1','G1','x1','A2','G2','x2','bg'},...
            'independent','x');
        opt=fitoptions(myfit);
        
        % Background is max
        bg=max(Y);
        
        % Find center
        [Ymin,ind]=min(Y);
        A=bg-Ymin;        
        xC=X(ind);
        
        % Assign guess
        G=[A 30 -160 A 30 -135 bg];
        
        
        opt.StartPoint=G;
        opt.Robust='bisquare';
        opt.Lower=[0 0 -inf 0 0 -inf 0];
        
        % Perform the fit
        fout=fit(X,Y,myfit,opt);
        disp(fout);

        % Plot the fit
        tt=linspace(min(X),max(X),1000);
        pF=plot(tt,feval(fout,tt),'r-','linewidth',1);
        lStr=['xC=(' num2str(round(fout.x1,1)) ',' num2str(round(fout.x2,1)) ')' ...
            ' FWHM=(' num2str(round(fout.G1,1)) ',' num2str(round(fout.G2,1)) ')' ];
        legend(pF,lStr,'location','best');
    end
    
    %% Triple Lorentzian
    if length(X)>4 && lorentz_triple
        myfit=fittype('bg+A1*(G1/2).^2*((x-x1).^2+(G1/2).^2).^(-1)+A2*(G2/2).^2*((x-x2).^2+(G2/2).^2).^(-1)+A3*(G3/2).^2*((x-x3).^2+(G3/2).^2).^(-1)',...
            'coefficients',{'A1','G1','x1','A2','G2','x2','A3','G3','x3','bg'},...
            'independent','x');
        opt=fitoptions(myfit);
        
        % Background is max
        bg=max(Y);
        
        % Find center
        [Ymin,ind]=min(Y);
        A=bg-Ymin;        
        xC=X(ind);
        
        % Assign guess
        G=[0.7 100 -220 0.7 100 -100 0.75 30 15 bg];
        
        
        opt.StartPoint=G;
        opt.Robust='bisquare';
        opt.Lower=[0 0 -inf 0 0 -inf 0];
        
        % Perform the fit
        fout=fit(X,Y,myfit,opt);
        disp(fout);

        % Plot the fit
        tt=linspace(min(X),max(X),1000);
        pF=plot(tt,feval(fout,tt),'r-','linewidth',1);
        lStr=['xC=(' num2str(round(fout.x1,1)) ',' num2str(round(fout.x2,1)) ',' num2str(round(fout.x3,1)) ')' ...
            ' FWHM=(' num2str(round(fout.G1,1)) ',' num2str(round(fout.G2,1)) ',' num2str(round(fout.G3,1)) ')' ];
        legend(pF,lStr,'location','best');
    end
    
    %% Four Lorentzian
    
    % Assymetric lorentzian fit, good for AM spec
    if length(X)>10 && fit_lorentz_assymetric_4
        g=@(x,a,x0,G) 2*G./(1+exp(a*(x-x0)));
        y=@(x,a,x0,G,A) A./(4*(x-x0).^2./g(x,a,x0,G).^2+1);    
        
        yTot = @(pa1,pa2,pa3,pa4,...
                pb1,pb2,pb3,pb4,...
                pc1,pc2,pc3,pc4,...
                pd1,pd2,pd3,pd4,...
                bg,x) ...
            y(x,pa1,pa2,pa3,pa4) + ...
            y(x,pb1,pb2,pb3,pb4) + ...
            y(x,pc1,pc2,pc3,pc4) + ...
            y(x,pd1,pd2,pd3,pd4) + bg;           
        
        myfit=fittype(@(pa1,pa2,pa3,pa4,...
                pb1,pb2,pb3,pb4,...
                pc1,pc2,pc3,pc4,...
                pd1,pd2,pd3,pd4,...
                bg,x) yTot(pa1,pa2,pa3,pa4,...
                pb1,pb2,pb3,pb4,...
                pc1,pc2,pc3,pc4,...
                pd1,pd2,pd3,pd4,...
                bg,x),'coefficients',{'pa1','pa2','pa3','pa4',...
                'pb1','pb2','pb3','pb4',...
                'pc1','pc2','pc3','pc4',...
                'pd1','pd2','pd3','pd4','bg'},...
            'independent','x'); 
        
        opt=fitoptions(myfit);
        
        opt.StartPoint = zeros(1,17);
        % assymetry
        opt.StartPoint(1) = 0;
        opt.StartPoint(5) = 0;
        opt.StartPoint(9) = 0;
        opt.StartPoint(13) = 0;
        
        % center
        opt.StartPoint(2) = 0;
        opt.StartPoint(6) = 50;
        opt.StartPoint(10) = 120;
        opt.StartPoint(14) = 190;
        
        % linewidth
        opt.StartPoint(3) = 30;
        opt.StartPoint(7) = 30;
        opt.StartPoint(11) = 30;
        opt.StartPoint(15) = 30;

        % ampltiude
        opt.StartPoint(4) = -2e4;
        opt.StartPoint(8) = -1e4;
        opt.StartPoint(12) = -1e4;
        opt.StartPoint(16) = -1e4;
        
        % bkacground
        opt.StartPoint(17) = 4E4;        
        
        G0=30;
        bg=min(Y);max(Y);
        A1=(max(Y)-min(Y));
        inds=[Y>.9*max(Y)];            
        
        [~,i]=max(Y);
        x0=X(i);
        x0=mean(X(inds));     
%         opt.StartPoint=[.1 -110 G0 A0 bg];  
        opt.Robust='bisquare';
%         opts.Weights=w;
        
        fout_lorentz=fit(X,Y,myfit,opt);
%         ci = confint(fout_lorentz,0.95);   
        
        XF=linspace(min(X)-5,max(X)+5,1000);
        xlim([min(X)-0.1 max(X)+0.1]);
        pExp=plot(XF,feval(fout_lorentz,XF),'r-','linewidth',2);
%         str=['$f_0 = ' num2str(round(fout_lorentz.x0,2)) '\pm' num2str(round((ci(2,2)-ci(1,2))/2,2)) '$ kHz' newline ...
%             '$\mathrm{FWHM} = ' num2str(round(abs(fout_lorentz.G),2)) ' $ kHz'];
%         legend(pExp,{str},'interpreter','latex','location','best','fontsize',8);         
    end
    
    %% Assymetric Lorentzian
    
    % Assymetric lorentzian fit, good for AM spec
    if length(X)>4 && lorentz_asym_single
        g=@(x,a,x0,G) 2*G./(1+exp(a*(x-x0)));
        y=@(x,a,x0,G,A,bg) A./(4*(x-x0).^2./g(x,a,x0,G).^2+1);        
        
        myfit=fittype(@(bg,a1,x1,G1,A1,x) y(x,a1,x1,G1,A1)+bg,...
            'coefficients',{'bg','a1','x1','G1','A1'},...
            'independent','x'); 
        
        opt=fitoptions(myfit);
        
        % Background
        bg=min(Y);max(Y);
        
        % Linewidth
        G1=30;

        % Contrast
        A1=(max(Y)-min(Y));
        
        % Center Point
        inds=[Y>.9*max(Y)];         
        x1=mean(X(inds));    
        
        % Assymetry
        a1 = -0.05; % Long on right
%         a1 = +0.05; % Long on left
        
        opt.StartPoint=[bg a1 x1 G1 A1];  
        opt.Robust='bisquare';
        
        fout_lorentz=fit(X,Y,myfit,opt)
        ci = confint(fout_lorentz,0.95);   
        
        XF=linspace(min(X)-5,max(X)+5,1000);
        xlim([min(X)-0.1 max(X)+0.1]);
        pExp=plot(XF,feval(fout_lorentz,XF),'r-','linewidth',2);
        str=['$f_0 = ' num2str(round(fout_lorentz.x1,2)) '\pm' ...
            num2str(round((ci(2,2)-ci(1,2))/2,2)) '$ kHz' newline ...
            '$\mathrm{FWHM} = ' ...
            num2str(round(abs(fout_lorentz.G1),2)) ' $ kHz' newline ...
            '$a = ' num2str(round(fout_lorentz.a1,3)) '$ kHz'];
        legend(pExp,{str},'interpreter','latex','location','northwest','fontsize',10);         
    end
    
    % Assymetric lorentzian fit
    if length(X)>9 && lorentz_asym_double
        g=@(x,a,x0,G) 2*G./(1+exp(a*(x-x0)));
        y=@(x,a,x0,G,A) A./(4*(x-x0).^2./g(x,a,x0,G).^2+1);   
        
        foo = @(x,bg,a1,x1,G1,A1,a2,x2,G2,A2) ...
            y(x,a1,x2,G1,A1)+y(x,a2,x2,G2,A2)+bg;
        
        myfit=fittype(@(bg,a1,x1,G1,A1,a2,x2,G2,A2,x) ...
            y(x,a1,x1,G1,A1)+y(x,a2,x2,G2,A2)+bg,...
            'coefficients',{'bg',...
            'a1','x1','G1','A1',...
            'a2','x2','G2','A2'},...
            'independent','x'); 
        opt=fitoptions(myfit);

        bg=min(Y);
        
        G1 = 13;
        G2 = 13;        
                       
        A1 = (max(Y)-min(Y));   
        A2 = A1/40;
        
                
        inds=[Y>.9*max(Y)];
        x1 =0; mean(X(inds)); 
        x2 = 10; %x2 = -155;
        
         A2 = (max(Y)-min(Y));   
        A1 = A2/30;
        
                
        inds=[Y>.9*max(Y)];
        x1 =0; mean(X(inds)); 
        x2 = 15; %x2 = -155;

        
        a1 = -.1;
        a2 = -.05;
        
        
        opt.StartPoint=[bg a1 x1 G1 A1 a2 x2 G2 A2];  
        opt.Robust='bisquare';
%         opts.Weights=w;
        
        fout_lorentz=fit(X,Y,myfit,opt)
        ci = confint(fout_lorentz,0.95);   
        
        XF=linspace(min(X)-5,max(X)+5,1000);
%         xlim([min(X)-0.1 max(X)+0.1]);
        pFit=plot(XF,feval(fout_lorentz,XF),'r-','linewidth',2);
        str=['$f_1 = ' num2str(round(fout_lorentz.x1,2)) '\pm' num2str(round((ci(2,3)-ci(1,3))/2,2)) '$ kHz' newline ...
            '$\mathrm{FWHM} = ' num2str(round(abs(fout_lorentz.G1),2)) ' $ kHz' newline ...
            '$f_2 = ' num2str(round(fout_lorentz.x2,2)) '\pm' num2str(round((ci(2,7)-ci(1,7))/2,2)) '$ kHz' newline ...
            '$\mathrm{FWHM} = ' num2str(round(abs(fout_lorentz.G2),2)) ' $ kHz' newline...
            'A1 =' num2str(round(fout_lorentz.A1,2)) newline...
            'A2 =' num2str(round(fout_lorentz.A2,2))];
        legend(pFit,str,'location','best','interpreter','latex');
    end
    
    %% Lorentzian
    
    if length(X)>4 && lorentz_single
        % Symmetric Lorentzian
        myfit=fittype('A*(G/2).^2*((x-x0).^2+(G/2).^2).^(-1)+bg','coefficients',{'A','G','x0','bg'},...
            'independent','x');
        opt=fitoptions(myfit);
        G0=20;
        bg=min(Y);
        A1=(max(Y)-min(Y));
        inds=[Y>.8*max(Y)];
        x0=mean(X(inds));
        opt.StartPoint=[A1/100 G0 x0 bg];   
%         opt.Upper=[1 3*G0 x0+range(X) 0];   

        opt.Robust='bisquare';


        fout_lorentz=fit(X,Y,myfit,opt);

        XF=linspace(min(X),max(X),1000);
        pExp=plot(XF,feval(fout_lorentz,XF),'r-','linewidth',2);

        str=['$f_0 = ' num2str(round(fout_lorentz.x0,2)) '$ kHz' newline ...
            '$\mathrm{FWHM} = ' num2str(round(abs(fout_lorentz.G),2)) ' $ kHz'];
        legend(pExp,{str},'interpreter','latex','location','best','fontsize',8);        
%         xlim([130 200]);    
    end
    
    %% Gauss Single
    if length(X)>4 && gauss_single
        myfit=fittype('bg+A1*exp(-(x-x1).^2/G1.^2)',...
            'coefficients',{'A1','G1','x1','bg'},...
            'independent','x');
        opt=fitoptions(myfit);
        % Background is max
        bg=min(Y);
        % Find center
        [Ymin,ind]=min(Y);
        A=bg-Ymin;
        xC=X(ind);
 
        G=[A 10 0 bg];
        opt.StartPoint=G;
        opt.Robust='bisquare';
%         opt.Lower=[0 0 -inf 0 0 -inf 0];
        % Perform the fit
        fout=fit(X,Y,myfit,opt);
        disp(fout);
        ci = confint(fout,0.95);
        
        % Plot the fit
        tt=linspace(min(X),max(X),1000);
        pF=plot(tt,feval(fout,tt),'r-','linewidth',1);
        lStr=['xC=(' num2str(round(fout.x1,2)) '±' ...
            num2str(abs(round(ci(1,3)-fout.x1,2))) ','...
             ')' ...
            ' FWHM=(' num2str(round(fout.G1,1)) ')' ];
        legend(pF,lStr,'location','best');
    end
    
    %% Rabi
    
    if length(X)>4 && Rabi_oscillation       
        
        guess_freq = 1/.08;
        guess_tau = 0.5;
%     
%         myfunc=@(N0,f,tau,t) N0*(1 - exp(-pi*t/tau).*cos(2*pi*f*t))/2;           
%         fitFuncStr = '$0.5N_0\left(1-\exp(-\pi t / \tau)\cos(2 \pi f t)\right)$';

    myfunc=@(N0,f,tau,t) N0*(1 - exp(-pi*t/tau).*cos(2*pi*f*t+pi))/2;   
    fitFuncStr = '$0.5N_0\left(1-\exp(-\pi t / \tau)\cos(2 \pi f t)\right)$';


    % Define the fit
    myfit=fittype(@(N0,f,tau,t) myfunc(N0,f,tau,t),'independent','t',...
        'coefficients',{'N0','f','tau'});
    opt=fitoptions(myfit);   
    
    
    opt.StartPoint=[max(Y) guess_freq guess_tau];
    opt.Lower=[max(Y)/5 .1 0];
    opt.Upper=[max(Y) 100 1000];

    opt.Robust='bisquare';

    % Perform the fit
    fout=fit(X,Y,myfit,opt);
    % Construct fit strings
    omega_rabi=2*pi*fout.f; 
    paramStr=['$N_0=' num2str(fout.N0,2) ',~f=' num2str(round(fout.f,2)) ...
        '~\mathrm{kHz},~\tau=' num2str(round(fout.tau,2)) '~\mathrm{ms}' ...
        '$'];

    tt=linspace(0,max(X),1000);
     pF=plot(tt,feval(fout,tt),'r-','linewidth',1);
    
    text(.45,.90,fitFuncStr,'units','normalized','interpreter','latex',...
        'horizontalalignment','right','fontsize',14);
    
    xL=get(gca,'XLim');
    yL=get(gca,'YLim');

    xlim([0 xL(2)]);
    ylim([0 yL(2)+.1]);

    legend(pF,paramStr,'location','northeast','interpreter','latex');
    outdata.Fit=fout;

    end
   
%% Power law
    if length(X)>4 && UniPowerLaw

        myfit=fittype('A1*(1+ t/(6*tau))^(-A0)',...
            'coefficients',{'A0','A1','tau'},...
            'independent','t');

        % Fit options and guess
        opt=fitoptions(myfit);        
        Ag = max(Y);
        A0 = 6;
        taug = median(X);
        G=[A0 Ag taug];        
        opt.StartPoint=G;
        opt.Lower = [0 0 0];

        % Perform the fit
        fout=fit(X,Y,myfit,opt);

        % Plot the fit
        tt=linspace(0,max(X),1000);
        xlim([0 max(X)]);
        pF=plot(tt,feval(fout,tt),'r-','linewidth',1);
        lStr=['$ \tau = ' num2str(round(fout.tau,3)) '~\mathrm{ms}$' ...
            newline ...
            '$A_0 = ' num2str(round(fout.A0,3),'%.3e') '$' ... 
            newline ...
            '$A_1 = ' num2str(round(fout.A1,3),'%.3e') '$'];
        legend(pF,lStr,'location','best','interpreter','latex');        
        str = '$A_0+ A_1\exp(-t/\tau)$';
        t=text(.02,.03,str,'units','normalized',...
            'fontsize',10,'interpreter','latex');
    end
   
%% Save the Figure

if doSave
    saveFigure(hFB,figName,saveOpts);
end
if doSave
    save([saveDir filesep 'custom_data'],'custom_data');
end
    
