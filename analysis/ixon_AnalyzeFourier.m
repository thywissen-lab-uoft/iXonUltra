
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


%% Animate cloud FFT
% ixon_doAnimateFFT = 1;
% 
% ixon_animateOptsFFT=struct;   
% 
% 
% % Variable to animate versus
% ixon_animateOptsFFT.xUnit=ixon_unit;
% 
% % Animation Timings
% ixon_animateOptsFFT.StartDelay=2; % Time to hold on first picture
% ixon_animateOptsFFT.MidDelay=.25;     % Time to hold in middle picutres
% ixon_animateOptsFFT.EndDelay=2;     % Time to hold final picture
% 
% % Animate in ascending or descending order?
% % animateOpts.Order='descend';    
% ixon_animateOptsFFT.Order='ascend';
% 
% % Color limit for image
% ixon_animateOptsFFT.CLim=[0 .5];   % Color limits 
% ixon_animateOptsFFT.CLim=[0 3];   % Color limits 
% 
% ixon_animateOptsFFT.CLim='auto';   % Automatically choose CLIM?
% 
% % FFT UV Cutoff
% % Reduce the animation view to within a frequency of 1/L
% ixon_animateOptsFFT.mask_UV=1;
% ixon_animateOptsFFT.LMin=20;
% 
% % FFT IR Cutoff
% % Apply mask to interior regions to mask 
% ixon_animateOptsFFT.mask_IR=1;
% ixon_animateOptsFFT.LMax=200;
% 
% if ixon_doAnimateFFT == 1 && ixon_doFFT && ixon_doSave
%     ixon_animateFFT(ixondata,ixon_xVar,ixon_animateOptsFFT);
% end