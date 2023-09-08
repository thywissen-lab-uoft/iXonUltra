
%% ANALYSIS : FFT BOX COUNT

% if ixon_fft_doBoxCount && ixon_doFFT    
%     fft_boxOpts=struct;
%     fft_boxOpts.maskIR=fft_opts.maskIR;
%     fft_boxOpts.LMax=fft_opts.LMax;
%     fft_boxOpts.maskUV=fft_opts.maskUV;
%     fft_boxOpts.LMin=fft_opts.LMin;
%         ixondata=ixon_fft_boxCount(ixondata,fft_boxOpts);
% end

%% PLOTTING : FFT BOX COUNT
% ixon_fft_boxPopts=struct;
% ixon_fft_boxPopts.xUnit=ixon_unit;
% if ixon_fft_doBoxCount  && ixon_doFFT 
%     % Plot the second moments
%     hF_ixon_size=ixon_fft_showBoxMoments(ixondata,ixon_xVar,ixon_fft_boxPopts);   
%     if ixon_doSave;ixon_saveFigure(ixondata,hF_ixon_size,'ixon_fft_box_size');end          
% end