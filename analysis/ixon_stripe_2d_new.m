%% 2D Stripe Analysis

% DataStuff

opt = struct;
opt.fig = figure;

for kk=1:length(ixondata)
   z=ixondata(kk).Z;
   stripeFit2D(z,opt) 
   
end