%% 2D Stripe Analysis
% This analzyes the stripes for take from the ixon camera.  This is a 2D
% fit over the entire cloud.  This fits a 2D gaussian modulated by a sine
% wave at a particular angle.  This is pariticularly useful to fit the
% angular dependence of the data. 
%
% The input fit parameters are specified in the options structure.


stripe_2d_opts=struct;

stripe_2d_opts.xUnit=ixon_unit;

stripe_2d_opts.ShimFit=0;
stripe_2d_opts.Theta=[-90 90]; % Specify the domain (MUST BE 180 DEGREES)
stripe_2d_opts.saveAnimation=1;        % save the animation?
stripe_2d_opts.StartDelay=1;
stripe_2d_opts.MidDelay=.5;
stripe_2d_opts.EndDelay=1;

if do_2dStripeAnalysis
    [hF_stripe_2d,stripe_data2d]=analyzeStripes2(ixondata,ixon_xVar,stripe_2d_opts);

    if ixon_doSave
        ixon_saveFigure(ixondata,hF_stripe_2d,'ixon_field_stripe_2d');        
    end
    
    outdata.stripe_data2d=stripe_data2d;   
end