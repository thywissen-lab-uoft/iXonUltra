  % Do basic analysis on raw counts
    ixondata=ixon_computeRawCounts(ixondata);

    % Plot histogram of raw counts
    hist_opts=struct;    
    hist_opts.xUnit=ixon_unit;

    % Specify the plot variable and units    
    hist_opts.Outliers=[10 50]; % Histogram wont plot outliers of this many low/high
    hist_opts.GlobalLimits=1;   % Maintain historgram x limits
    hist_opts.BinWidth=10;       % Histogram bin width
    hist_opts.ImageNumber=1;    % Which image to histogram (overwritten)
    hist_opts.YScale='Log';     % Histogram y scale
    % hist_opts.YScale='Linear';

    for kk=1:size(ixondata(1).RawImages,3)
        hist_opts.ImageNumber=kk;
        hF_ixon_rawhist=ixon_showRawCountHistogram(ixondata,ixon_xVar,hist_opts);
        if ixon_doSave;ixon_saveFigure(ixondata,hF_ixon_rawhist,['ixon_raw_hist' num2str(kk)]);end
    end
    % Plot raw count total
    raw_opts=struct;   
    raw_opts.xUnit=ixon_unit;    
    % Define the variable and units    
    raw_opts.FitLinear=0;    
    hF_ixon_rawtotal=ixon_showRawCountTotal(ixondata,ixon_xVar,raw_opts);
    if ixon_doSave;ixon_saveFigure(ixondata,hF_ixon_rawtotal,['ixon_raw_counts']);end