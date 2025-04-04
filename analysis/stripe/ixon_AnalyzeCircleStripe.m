function ixondata = ixon_AnalyzeCircleStripe(ixondata)
    opts = struct;
    opts.doDebug = 0;
    for kk=1:length(ixondata)
        disp([num2str(kk) '/' num2str(length(ixondata)) ' circle stripe'])
        X = ixondata(kk).X;
        Y = ixondata(kk).Y;
        Z = ixondata(kk).Z;
        myz=sum(Z,3);    
        stripe = StripeCircle(X,Y,myz,opts);
        stripe.Params= ixondata(kk).Params;        
        ixondata(kk).StripeCircular = stripe;
    end 
end

