function ixondata = ixon_fft_maskIR(ixondata,LMax)
%% Add IR Cuttoff

% Make the mask
F_IR=1/LMax;
X=ixondata(1).fft_F;
Y=ixondata(1).fft_F;
[xx,yy]=meshgrid(X,Y);
rr=sqrt(xx.^2+yy.^2);
IR_mask=rr>F_IR;

for kk=1:length(ixondata)
   ixondata(kk).fft_Z= ixondata(kk).fft_Z.*IR_mask;
end
 

end

