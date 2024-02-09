% ixon_main_qgm.m
%
% Author : CF Fujiwara
%
% This script is the primary analysis code for ixon images which are single
% plane.

% Display this filename
disp(repmat('-',1,60));disp(repmat('-',1,60));    
disp(['Calling ' mfilename '.m']);
disp(repmat('-',1,60));disp(repmat('-',1,60));  

%% Bin Stripe
if ixon_doQGM_BinStripe
    LGuess = 25;
    ColorThreshold = [1000 3000];
    
    if ~isfield(qgmdata,'LatticeBin')
        return;
    end
    clear out
    for n = 1:length(qgmdata)                
        for kk = 1:length(qgmdata(n).LatticeBin)
            n1 = qgmdata(n).LatticeBin(kk).n1;
            n2 = qgmdata(n).LatticeBin(kk).n2;
            Zb = qgmdata(n).LatticeBin(kk).Zbin;    
            opts_stripe.LGuess = LGuess;
            opts_stripe.FigNum=3000+10*(n-1)+kk-1;
            opts_stripe.FigNum=3000;

            opts_stripe.ColorThreshold = ColorThreshold;
            
            [out(kk),hF_bin_stripe] = ixon_fitStripe_dig(n1,n2,Zb,opts_stripe);
%                   exportgraphics(gcf,'testAnimated.gif','Append',true);
            frame=getframe(hF_bin_stripe);
            im = frame2im(frame);
            [A,map] = rgb2ind(im,256);  
            
            filename = fullfile(saveDir,'binstripe_75.gif');
            
            if out(kk).ModDepth>=.75 && out(kk).Counts>0.5e6
                switch n
                    case 1
                        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
                    case length(qgmdata)
                        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1);
                    otherwise
                        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',.1);
                end
            end

        end
        qgmdata(n).BinStripe = out;    
        qgmdata(n).BinStripe = out;
    end
end  

%% Digitization Stuff
if ixon_doQGM_Digitize
    ixon_main_digital;
end