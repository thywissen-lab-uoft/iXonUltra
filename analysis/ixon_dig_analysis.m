% ixon_main_dig.m
%
% Author : CF Fujiwara
%
% This script is the primary analysis code for ixon images which are single
% plane.

% Display this filename
disp(repmat('-',1,60));disp(repmat('-',1,60));    
disp(['Calling ' mfilename '.m']);
disp(repmat('-',1,60));disp(repmat('-',1,60));    

% Add all subdirectories for this m file
curpath = fileparts(mfilename('fullpath'));
addpath(curpath);addpath(genpath(curpath));

a = fileparts(curpath);
addpath(a);addpath(genpath(a));
%% Close all non GUI figures
% Close all figures without the GUI tag.
figs=get(groot,'Children');
disp('Closing all non GUI figures.');
for kk=1:length(figs)
   if ~isequal(figs(kk).Tag,'GUI')
        disp(['Closing figure ' num2str(figs(kk).Number) ' ' figs(kk).Name]);
        close(figs(kk)) 
   end
end
disp(' ');
%% Select image directory
    
if ~exist('dig_auto_file')
   dig_auto_file = 1; 
end

if dig_auto_file
    dialog_title='Select GUI data';       
    [filename,dig_imgdir,b]=uigetfile(fullfile(ixon_getDayDir,'*.mat'),dialog_title);
    filename = fullfile(dig_imgdir,filename);    
    if  b == 0
        disp('Canceling.');    
        return; 
    end
end
%% Load the Data
clear digdata
tic
fprintf('Loading digdata ...');
digdata = load(filename);
disp([' done (' num2str(toc,'%.2f') 's)']);

%% Initialize Options

dig_opts = struct;
dig_opts.Quality = 'auto';
dig_opts.saveDir=dig_imgdir;
dig_opts.FigLabel=digdata.SourceDirectory{1};

%% Analysis Variable
% This section of code chooses the variable to plot against for aggregate
% plots.  The chosen variable MUST match a variable provided in the params
% field of the .mat file. The unit has no tangibile affect and only affects
% display properties.

% Choose what kind of variable to plot against (sequencer/camera)
dig_opts.varType        = 'param';          % always select 'param' for now 
dig_opts.autoXVar       = 0;                % Auto detect changing variable?
dig_opts.autoUnit       = 1;                % Auto detect unit for variable?
dig_opts.xVar           = 'conductivity_mod_time';  % Variable Name
dig_opts.xVar           = 'ExecutionDate';  % Variable Name
dig_opts.overrideUnit   = 'V';              % If ixon_autoUnit=0, use this
dig_opts.doSave         = 1;                % Save Analysis?

%% Flags

% Recenter all binned data to have same limits
dig_doShowCloud                         = 1;
dig_doShowCloudAnimate                  = 1;
dig_standardAnalysis                    = 1;
dig_ac_conductivity_fit                 = 0;
dig_doRadialAnalysis                    = 1;
%% Show CLoud

if dig_doShowCloud
    opts = dig_opts;
    opts.doAnimate = 1;
    opts.ROI = 'max';
    opts.filename = 'dig_animateFull.gif';
    hF_digCloud = dig_showCloud(digdata,opts);

    if dig_opts.doSave
        ixon_saveFigure2(hF_digCloud,...
         'dig_average',dig_opts);  
    end
end

%% Standard Analysis

if dig_standardAnalysis
    hF_digStandard = dig_showStandardAnalysis(digdata,opts);
    if dig_opts.doSave
        ixon_saveFigure2(hF_digStandard,...
         'dig_standard',dig_opts);  
    end
end

%% Radial Analysis

if dig_doRadialAnalysis
    opts = dig_opts;   
    opts.BinStep = 3;                   % delta r bin size    

    opts.useAverageCenter = 0;
    digdata = dig_compute_radial(digdata,opts);     % Compute radial profile for all images

    
    opts.rMaxShow = 80;                 % max r to plot
    opts.nMaxShow = 0.5;               % max density to plot
    opts.showDevParametrization  = 0;   % show standard deviation?
    
    % Show radial profiles
    hFs=dig_showRadialProfile(digdata,opts);     
     if dig_opts.doSave
         for kk=1:length(hFs)
             ixon_saveFigure2(hFs(kk),...
                ['dig_radial_profile_' num2str(kk)],dig_opts);  
         end
     end        
     
      opts.ForceAverage = 0;
     if isequal(digdata.xVar,'ExecutionDate')
        opts.ForceAverage = 1; 
     end
%     opts.ForceAverage = 0;
    opts.doAnimate = 1;
    
    [hFs_radial] = dig_radialAnalysis_average_images(digdata,opts);
    if dig_opts.doSave
       for kk=1:length(hFs_radial)
           ixon_saveFigure2(hFs_radial(kk),...
                ['dig_radial_' num2str(kk)],dig_opts);   
       end
    end
    
    doExportforDrut=0;
    if doExportforDrut    
       % Constants
        h = 6.62607015*10^-34;
        m = 39.96399848*1.66053906660*10^-27;
        lam = (1054*10^-9*1054*10^-9*1064*10^-9)^(1/3);
        aL = lam/2;

        % 120mW XDT + 2.5ER request lattice trap frequencies - calibrated
        % 02/29/24 (reanalyzed 04/11/24)
        omega_x = 2*pi*67.3;
        omega_y = 2*pi*60.2;
        omega_bar = (omega_x*omega_y).^(1/2);

        % Harmonic potential in Hz
        PotentialVector = 0.5*m*omega_bar^2*aL^2.*rVec.^2/h;

        kappa = 0.5*m*omega_bar^2*aL^2/h;

        % Assign to output
        out = struct;
        out.Images                         = digdata.Z;
        out.SiteVector1                    = digdata.n1;
        out.SiteVector2                    = digdata.n2;            
        out.RadialVector                   = digdata.r;
        out.RadialOccupation               = digdata.Zr;
        out.RadialPotential                = digdata.r.^2*kappa;                  

        try if ~exist(dig_opts.saveDir,'dir');mkdir(dig_opts.saveDir);end;end
        filename = fullfile(dig_opts.saveDir,'dig_radial_data.mat');
        disp(['Saving ' filename ' ...']);
        save(filename, '-struct','dig_radial_data');
    end
     
    % obsolte analysis
%     [hF_digRadial_2,dig_radial_data] = dig_radialAnalysis(digdata,opts);        


end

%% Conductivity Analysis

if dig_ac_conductivity_fit
opts = dig_opts;

[hF_conductivity,conductivity_data] = dig_ac_conductivity(digdata,opts);
    if dig_opts.doSave
        ixon_saveFigure2(hF_conductivity,...
         'dig_conductivity',dig_opts);  
    end
    try if ~exist(dig_opts.saveDir,'dir');mkdir(dig_opts.saveDir);end;end
    filename = fullfile(dig_opts.saveDir,'conductivity_data.mat');
    disp(['Saving ' filename ' ...']);
    save(filename, '-struct','conductivity_data');


end

% %%
% Xc2vec = zeros(1,size(digdata.Zdig,3));
% 
% 
% Nmed = median(sum(digdata.Zdig,[1 2]));
% % 
% % for nn=1:size(digdata.Zdig,3)
% %     y=smooth(sum(digdata.Zdig(:,:,nn),2),10);
% %     
% %     try
% %     [pks,locs,w,p]=findpeaks(y,'SortStr','descend');
% %     
% %     pks=pks(1:3);
% %     locs=locs(1:3);
% %     [val,ind]=min(abs(locs-75));
% %     loc = locs(ind);
% %     
% %     n2sub_ind = loc + [-12:1:12];
% %     
% %     Zsub = digdata.Zdig(n2sub_ind,:,nn);    
% %     n1sub = digdata.n1;
% %     n2sub = digdata.n2(n2sub_ind);
% %     
% %     [n1,n2] = meshgrid(n1sub,n2sub);
% %     
% %     Xc2= sum(Zsub.*n1,'all')/sum(Zsub,'all');  
% %     catch 
% %         Xc2 = NaN;
% %     end
% %     Xc2vec(nn) = Xc2;    
% % end
% 
% inds=[sum(digdata.Zdig,[1 2])<Nmed/3];
% 
% % digdata.Xc75 =Xc2vec;
% 
% x=digdata.X'+150;
% % y = digdata.Xc75';
% y = digdata.Xc';
% % y = (digdata.Xc-(rx0+65)/.527)';
% % inds = isnan(y);
% 
% x(inds)=[];
% y(inds)=[];
% 
% 
% hF=figure(7002);
% clf
% hF.Color='w';
% hF.Position=[100 100 600 300];
% co=get(gca,'colororder');
% plot(x,y,'o','markerfacecolor',co(1,:),...
%     'markersize',8,'markeredgecolor',co(1,:)*.5);
% xlabel([digdata.xVar ' + 150 ms'],'interpreter','none');
% ylabel('x center (sites)');
% hold on
% 
% P=[digdata.Params];
% f = unique([P.conductivity_mod_freq]);
% myfit = fittype(@(A,phi,x0,t) A*sin(2*pi*f*t*1e-3+phi) + x0 ,...
%     'independent','t','coefficients',{'A','phi','x0'});
% % myfit = fittype(@(A,phi,x0,v0,a0,t) A*sin(2*pi*f*t*1e-3+phi) + x0 +v0*t +0.5*a0^2*t.^2 ,...
% %     'independent','t','coefficients',{'A','phi','x0','v0','a0'});
% Ag = max(y)-min(y);
% Bg = Ag;
% % x0g = median(digdata.Xc75);
% fitopt = fitoptions(myfit);
% fitopt.StartPoint = [0.5*Ag 3.5 mean(y);];
% 
% % fitopt.StartPoint = [0.5*Ag 3.5 mean(y) 0 0];
% % fitopt.Upper = [10 10 inf 2 .02];
% % fitopt.Lower = [0 0 -inf -2 -.02];
% 
% fitopt.Robust = 'bisquare';
% % fitopt.Lower = [0 -20 x0g-5];
% 
% 
% fout = fit(x,y,myfit,fitopt)
% 
% str = ['A: ' num2str(round(fout.A,4)) ];
% P = [digdata.Params];
% B = unique([P.conductivity_FB_field]);
% Pevap = unique([P.xdt_evap1_power]);
% amp = unique([P.conductivity_ODT2_mod_amp]);
% 
% str = [num2str(B) ' G, ' num2str(Pevap*1e3) ' mW, ' ...
%     num2str(f) ' Hz, ' num2str(amp) ' V amp'];
% 
% title(str);
% % 
% cc=confint(fout,0.66);
% 
% Aerr =  (cc(2,1)-cc(1,1))/2;
% phierr = (cc(2,2)-cc(1,2))/2;
% x0err = (cc(2,3)-cc(1,3))/2;
% 
% str2 = ['A=' num2str(round(fout.A,2)) '\pm' num2str(round(Aerr,2)) ...
%     ' \phi=' num2str(round(fout.phi,2)) '\pm' num2str(round(phierr,2))];
% 
% % tt=linspace(200,200+1e3*3/f,1e3);
% % f=55;
% % phi=2.5;
% % x0=106.3;
% % A=2;
% % plot(tt,A*sin(2*pi*f*tt*1e-3+phi)+x0,'k-');
% plot(fout);
% 
% legend(str2,'location','best');
% s3 = digdata.SourceDirectory{1};
% t=uicontrol('style','text','String',s3,'backgroundcolor','w',...
%     'fontsize',6);
% t.Position(3:4)=t.Extent(3:4);
% t.Position(1:2)=[0 0];
% xlabel('time (ms)');
% 
% if dig_opts.doSave
%         ixon_saveFigure2(hF,...
%          'dig_center_stripe_fit',dig_opts);  
%     end