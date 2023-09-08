%% Fitting Options and Plotting Options
plt_opts = struct;
plt_opts.xUnit = ixon_unit;
plt_opts.xVar = ixon_xVar;
plt_opts.NumberScale = 'Linear';
% plt_opts.NumberScale = 'Log';
plt_opts.PositionUnit = 'px';
% plt_opts.PositionUnit = 'um';

% Number of counts fit
plt_opts.FitNumber_Exp          = 0;
plt_opts.FitNumber_Exp2Sum      = 0;
plt_opts.FitNumber_Lorentzian   = 0;

% Center position fit
plt_opts.FitCenter_Sine         = 0;
plt_opts.FitCenter_SineDecay    = 0;
plt_opts.FitCenter_SineGrow     = 0;
plt_opts.FitCenter_Linear       = 0;
%% PLOTTING : BOX COUNT

ixon_boxPopts=struct;
ixon_boxPopts.xUnit=ixon_unit;
ixon_boxPopts.NumberScale='Linear';
ixon_boxPopts.NumberExpFit              = 0;
ixon_boxPopts.NumberExp2SumFit          = 0;
ixon_boxPopts.NumberLorentzianFit       = 0;
ixon_boxPopts.CenterSineFit             = 0;       % Fit sine fit to cloud center
ixon_boxPopts.CenterDecaySineFit        = 1;  % Fit decaying sine to cloud center
ixon_boxPopts.CenterGrowSineFit         = 0;  % Fit decaying sine to cloud center
ixon_boxPopts.CenterLinearFit           = 0;     % Linear fit to cloud center

if ixon_doBoxCount  
    % Plot the atom number
    [hF_ixon_numberbox,Ndatabox]=ixon_showBoxNumber(ixondata,ixon_xVar,ixon_boxPopts);      
    yl=get(gca,'YLim');
    set(gca,'YLim',[0 yl(2)]);
    if ixon_doSave;ixon_saveFigure(ixondata,hF_ixon_numberbox,'ixon_box_number');end     
    
    % Plot the second moments
    hF_ixon_size=ixon_showBoxMoments(ixondata,ixon_xVar,ixon_boxPopts);   
    if ixon_doSave;ixon_saveFigure(ixondata,hF_ixon_size,'ixon_box_size');end     
    
    % Plot the cloud center
    [hF_ixon_center,Xc,Yc]=ixon_showBoxCentre(ixondata,ixon_xVar,ixon_boxPopts); 
    if ixon_doSave;ixon_saveFigure(ixondata,hF_ixon_center,'ixon_box_centre');end    
end


%% PLOTTING : GAUSSIAN
ixon_gauss_opts.xUnit=ixon_unit;
ixon_gauss_opts.NumberExpFit = 0;        % Fit exponential decay to atom number
ixon_gauss_opts.NumberLorentzianFit=0;   % Fit atom number to lorentzian
ixon_gauss_opts.NumberScale = 'linear'; 
% ixon_gauss_opts.NumberScale = 'log'; 

ixon_gauss_opts.CenterSineFit = 0;       % Fit sine fit to cloud center
ixon_gauss_opts.CenterDecaySineFit = 0;  % Fit decaying sine to cloud center
ixon_gauss_opts.CenterParabolaFit = 0;
ixon_gauss_opts.CenterLinearFit = 0;     % Linear fit to cloud center

if ixon_doGaussFit
    % Statistics if no variable is changing
    if isequal(ixon_xVar,'ExecutionDate')
        hF_stats=ixon_showGaussStats(ixondata);     
        if ixon_doSave;ixon_saveFigure(ixondata,hF_stats,'ixon_gauss_stats');end
    end
       
    % Counts
    [hF_numbergauss,Ndatagauss]=ixon_showGaussNumber(ixondata,ixon_xVar,ixon_gauss_opts);  
%  xlim([0 1.4]);
    if ixon_doSave;ixon_saveFigure(ixondata,hF_numbergauss,'ixon_gauss_number');end    
    
    % Size
    hF_size=ixon_showGaussSize(ixondata,ixon_xVar,ixon_gauss_opts);    
    if ixon_doSave;ixon_saveFigure(ixondata,hF_size,'ixon_gauss_size');end
        
    % Aspect Ratio
    hF_ratio=ixon_showGaussAspectRatio(ixondata,ixon_xVar,ixon_gauss_opts);    
    if ixon_doSave;ixon_saveFigure(ixondata,hF_ratio,'ixon_gauss_ratio');end
    
    % Centre
    hF_Centre=ixon_showGaussCentre(ixondata,ixon_xVar,ixon_gauss_opts);    
    if ixon_doSave;ixon_saveFigure(ixondata,hF_Centre,'ixon_gauss_position');end
    
%     hF_stats=ixon_showGaussStats(ixondata,ixon_gauss_opts);     
%     if ixon_doSave;ixon_saveFigure(ixondata,hF_stats,'ixon_stats');end

        
     % Style of profile --> cut or sum?
    style='cut';
%     style='sum';
    clear hF_X;    
    clear hF_Y;
    hF_X=[];
    hF_Y=[];
    
    hF_Xs=ixon_showGaussProfile(ixondata,'X',style,ixon_xVar,ixon_gauss_opts);        
    hF_Ys=ixon_showGaussProfile(ixondata,'Y',style,ixon_xVar,ixon_gauss_opts);  

%   Save the figures (this can be slow)
    if ixon_doSave
        for kk=1:length(hF_Xs)            
            ixon_saveFigure(ixondata,hF_Xs(kk),['ixon_gauss_profile_X'  num2str(kk)]);
        end
        for kk=1:length(hF_Ys)
            ixon_saveFigure(ixondata,hF_Ys(kk),['ixon_gauss_profile_Y' '_' num2str(kk)]);
        end
    end    
end
