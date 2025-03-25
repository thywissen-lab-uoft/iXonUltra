%% 1D STRIPE ANALYSIS
% This analyzes the stripes for field stability.  This is a 1D analysis,
% where a 1D sine wave with a gaussain envelope is fitted over a sub ROI of
% the data.
%
% This one dimensional analysis requires a variety of inputs to fit.  See
% the "1D Fit Parameters" below.

stripe_opts=struct;

% X plot unit
strope_opts.xUnit=ixon_unit;

% 1D Fit Parameters
stripe_opts.theta=57;               % Rotation Angle
stripe_opts.rotrange=[220 300];     % Sub region to inspect
stripe_opts.FitType='Sine';         % Fit Type
stripe_opts.LowThreshold=0.2;       % Low ampltude to ignore
stripe_opts.L0=80;                 % Guess wavelength (pixels)
stripe_opts.phi0=pi/2;             % Guess phase (radians)
stripe_opts.B0=0.4;                 % Guess modulation

% Animate
stripe_opts.saveAnimation=1;        % save the animation?
stripe_opts.StartDelay=.25;
stripe_opts.MidDelay=.25;
stripe_opts.EndDelay=.25;

% Field Analysis
field_opts=struct;
field_opts.xUnit=ixon_unit;
field_opts.FieldGradient=210;   % In G/cm
field_opts.LatticeSpacing=532E-9; % in meter
field_opts.FitType='Exp';

if doStripeAnalysis
    [hF_stripe,stripe_data]=analyzeStripes(ixondata,ixon_xVar,stripe_opts);
    
    if ixon_doSave;ixon_saveFigure(ixondata,hF_stripe,'ixon_stripe');end

    field_gradient=210; % field gradient in G/cm
    stripe_data.grabMagnetometer=1;
    stripe_data.Nsmooth=10;    
    
    [hF_field_stripe1,hF_field_stripe2,hF_field_sense]=fieldAnalysis(stripe_data,field_opts);
    
    if ixon_doSave
        ixon_saveFigure(ixondata,hF_field_stripe1,'ixon_field_stripe1');  
        ixon_saveFigure(ixondata,hF_field_stripe2,'ixon_field_stripe2');        

         if stripe_data.grabMagnetometer
            ixon_saveFigure(ixondata,hF_field_sense,'ixon_field_sense');
        end  
    end
    
    outdata.stripe_data=stripe_data;   
end

