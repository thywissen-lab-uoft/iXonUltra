%% Plotting Options
ixon_plt_opts = struct;
ixon_plt_opts.FigLabel = FigLabel;
ixon_plt_opts.xUnit = ixon_unit;
ixon_plt_opts.xVar = ixon_xVar;
ixon_plt_opts.NumberScale = 'Linear';
% plt_opts.NumberScale = 'Log';
ixon_plt_opts.PositionUnit = 'px';
% plt_opts.PositionUnit = 'um';

%% Fitting Options
ixon_fit_opts = struct;

% Number of counts fit
ixon_fit_opts.NumberExpFit          = 0;
ixon_fit_opts.NumberExpOffsetFit      = 0;
ixon_fit_opts.Number_Lorentzian   = 0;

% Center position fit
ixon_fit_opts.Center_Sine         = 0;
ixon_fit_opts.Center_SineDecay    = 0;
ixon_fit_opts.Center_SineGrow     = 0;
ixon_fit_opts.Center_Linear       = 1;


%% PLOTTING : BOX

if ixon_doBoxCount  
    % Counts
    hF_ixon_box_counts = ixon_showCounts(ixon_boxdata,ixon_xVar,ixon_plt_opts,ixon_fit_opts);
    if ixon_doSave;ixon_saveFigure2(hF_ixon_box_counts,'ixon_box_counts',saveOpts);end     
    
    if ~isequal(ixon_xVar,'ExecutionDate')
        ixon_plt_opts2 = ixon_plt_opts;
        ixon_plt_opts2.xVar='ExecutionDate';
    hF_ixon_box_counts_time = ixon_showCounts(ixon_boxdata,'ExecutionDate',ixon_plt_opts2,ixon_fit_opts);

        hF_ixon_box_counts_time.Position(2)=700;
        ylim([0 max(get(gca,'YLim'))]);    
        if ixon_doSave;ixon_saveFigure2(hF_ixon_box_counts_time,'ixon_box_counts_time',saveOpts);end
    end
    
    
    % Size
    hF_ixon_box_size = ixon_showSize(ixon_boxdata,ixon_xVar,ixon_plt_opts,ixon_fit_opts);
    if ixon_doSave;ixon_saveFigure2(hF_ixon_box_size,'ixon_box_size',saveOpts);end     
    
    % Center
    hF_ixon_box_centre = ixon_showCentre(ixon_boxdata,ixon_xVar,ixon_plt_opts,ixon_fit_opts);
    if ixon_doSave;ixon_saveFigure2(hF_ixon_box_centre,'ixon_box_centre',saveOpts);end     
end

%% PLOTTING : GAUSS

if ixon_doGaussFit  
    % Counts
    hF_ixon_gauss_counts = ixon_showCounts(ixon_gaussdata,ixon_xVar,ixon_plt_opts,ixon_fit_opts);
    if ixon_doSave;ixon_saveFigure2(hF_ixon_gauss_counts,'ixon_gauss_counts',saveOpts);end     
    
    % Size
    hF_ixon_gauss_size = ixon_showSize(ixon_gaussdata,ixon_xVar,ixon_plt_opts,ixon_fit_opts);
    if ixon_doSave;ixon_saveFigure2(hF_ixon_gauss_size,'ixon_gauss_size',saveOpts);end     
    
    % Center
    hF_ixon_gauss_centre = ixon_showCentre(ixon_gaussdata,ixon_xVar,ixon_plt_opts,ixon_fit_opts);
    if ixon_doSave;ixon_saveFigure2(hF_ixon_gauss_centre,'ixon_gauss_centre',saveOpts);end     
end

%% PLOTTING : GAUSSIAN
% 
% if ixon_doGaussFit
%     % Statistics if no variable is changing
%     if isequal(ixon_xVar,'ExecutionDate')
%         hF_stats=ixon_showGaussStats(ixondata);     
%         if ixon_doSave;ixon_saveFigure(ixondata,hF_stats,'ixon_gauss_stats');end
%     end
%        
%     % Counts
%     [hF_numbergauss,Ndatagauss]=ixon_showGaussNumber(ixondata,ixon_xVar,ixon_plt_opts);  
%     if ixon_doSave;ixon_saveFigure(ixondata,hF_numbergauss,'ixon_gauss_number');end    
%     
%     % Size
%     hF_size=ixon_showGaussSize(ixondata,ixon_xVar,ixon_plt_opts);    
%     if ixon_doSave;ixon_saveFigure(ixondata,hF_size,'ixon_gauss_size');end
%         
%     % Aspect Ratio
%     hF_ratio=ixon_showGaussAspectRatio(ixondata,ixon_xVar,ixon_plt_opts);    
%     if ixon_doSave;ixon_saveFigure(ixondata,hF_ratio,'ixon_gauss_ratio');end
%     
%     % Centre
%     hF_Centre=ixon_showGaussCentre(ixondata,ixon_xVar,ixon_plt_opts);    
%     if ixon_doSave;ixon_saveFigure(ixondata,hF_Centre,'ixon_gauss_position');end
%     
% %     hF_stats=ixon_showGaussStats(ixondata,ixon_plt_opts);     
% %     if ixon_doSave;ixon_saveFigure(ixondata,hF_stats,'ixon_stats');end
% 
%         
%      % Style of profile --> cut or sum?
%     style='cut';
% %     style='sum';
%     clear hF_X;    
%     clear hF_Y;
%     hF_X=[];
%     hF_Y=[];
%     
%     hF_Xs=ixon_showGaussProfile(ixondata,'X',style,ixon_xVar,ixon_plt_opts);        
%     hF_Ys=ixon_showGaussProfile(ixondata,'Y',style,ixon_xVar,ixon_plt_opts);  
% 
% %   Save the figures (this can be slow)
%     if ixon_doSave
%         for kk=1:length(hF_Xs)            
%             ixon_saveFigure(ixondata,hF_Xs(kk),['ixon_gauss_profile_X'  num2str(kk)]);
%         end
%         for kk=1:length(hF_Ys)
%             ixon_saveFigure(ixondata,hF_Ys(kk),['ixon_gauss_profile_Y' '_' num2str(kk)]);
%         end
%     end    
% end
