%% Fitting Options and Plotting Options
ixon_plt_opts = struct;
ixon_plt_opts.xUnit = ixon_unit;
ixon_plt_opts.xVar = ixon_xVar;
ixon_plt_opts.NumberScale = 'Linear';
% plt_opts.NumberScale = 'Log';
ixon_plt_opts.PositionUnit = 'px';
% plt_opts.PositionUnit = 'um';

% Number of counts fit
ixon_plt_opts.FitNumber_Exp          = 0;
ixon_plt_opts.FitNumber_Exp2Sum      = 0;
ixon_plt_opts.FitNumber_Lorentzian   = 0;

% Center position fit
ixon_plt_opts.FitCenter_Sine         = 0;
ixon_plt_opts.FitCenter_SineDecay    = 0;
ixon_plt_opts.FitCenter_SineGrow     = 0;
ixon_plt_opts.FitCenter_Linear       = 0;


%% PLOTTING : BOX COUNT

if ixon_doBoxCount  
    % Plot the atom number
%     [hF_ixon_numberbox,Ndatabox]=ixon_showBoxNumber(ixondata,ixon_xVar,ixon_plt_opts);      
%     yl=get(gca,'YLim');
%     set(gca,'YLim',[0 yl(2)]);
%     if ixon_doSave;ixon_saveFigure(ixondata,hF_ixon_numberbox,'ixon_box_number');end     

    hF_counts = showCounts(ixon_boxdata,ixon_xVar,ixon_plt_opts);
    
    % Plot the second moments
    hF_ixon_size=ixon_showBoxMoments(ixondata,ixon_xVar,ixon_plt_opts);   
    if ixon_doSave;ixon_saveFigure(ixondata,hF_ixon_size,'ixon_box_size');end     
    
    % Plot the cloud center
    [hF_ixon_center,Xc,Yc]=ixon_showBoxCentre(ixondata,ixon_xVar,ixon_plt_opts); 
    if ixon_doSave;ixon_saveFigure(ixondata,hF_ixon_center,'ixon_box_centre');end    
end


%% PLOTTING : GAUSSIAN

if ixon_doGaussFit
    % Statistics if no variable is changing
    if isequal(ixon_xVar,'ExecutionDate')
        hF_stats=ixon_showGaussStats(ixondata);     
        if ixon_doSave;ixon_saveFigure(ixondata,hF_stats,'ixon_gauss_stats');end
    end
       
    % Counts
    [hF_numbergauss,Ndatagauss]=ixon_showGaussNumber(ixondata,ixon_xVar,ixon_plt_opts);  
%  xlim([0 1.4]);
    if ixon_doSave;ixon_saveFigure(ixondata,hF_numbergauss,'ixon_gauss_number');end    
    
    % Size
    hF_size=ixon_showGaussSize(ixondata,ixon_xVar,ixon_plt_opts);    
    if ixon_doSave;ixon_saveFigure(ixondata,hF_size,'ixon_gauss_size');end
        
    % Aspect Ratio
    hF_ratio=ixon_showGaussAspectRatio(ixondata,ixon_xVar,ixon_plt_opts);    
    if ixon_doSave;ixon_saveFigure(ixondata,hF_ratio,'ixon_gauss_ratio');end
    
    % Centre
    hF_Centre=ixon_showGaussCentre(ixondata,ixon_xVar,ixon_plt_opts);    
    if ixon_doSave;ixon_saveFigure(ixondata,hF_Centre,'ixon_gauss_position');end
    
%     hF_stats=ixon_showGaussStats(ixondata,ixon_plt_opts);     
%     if ixon_doSave;ixon_saveFigure(ixondata,hF_stats,'ixon_stats');end

        
     % Style of profile --> cut or sum?
    style='cut';
%     style='sum';
    clear hF_X;    
    clear hF_Y;
    hF_X=[];
    hF_Y=[];
    
    hF_Xs=ixon_showGaussProfile(ixondata,'X',style,ixon_xVar,ixon_plt_opts);        
    hF_Ys=ixon_showGaussProfile(ixondata,'Y',style,ixon_xVar,ixon_plt_opts);  

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
