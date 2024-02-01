
function Y=unwrapPhaseTime(X,Y,Navg)
% CF: I think this will be buggy if you have repititions whose noise is
% larger than a phase, but if that happens you're data is too noisy anyway.

% Unwrap the phase with a given X and Y.  Here we shall try to minimize the
% change in slope of the 

if nargin ==2
   Navg = 4; 
end

[X,inds] = sort(X);
Y = Y(inds);

% Performa a basic unwrap to keep jump less than 2*pi
Y = unwrap(Y);


Y=mod(Y,2*pi);
Y0 = Y(1);


foo = @(Y,Y0) (Y/(2*pi)-round((Y-Y0)/(2*pi)))*2*pi;

Y=foo(Y,Y0);


% Navg = 4;
for kk=(1+Navg):length(Y)
%     keyboard
    Y0 = median(Y(kk-Navg:kk-1));
    Y(kk) = foo(Y(kk),Y0);
    
%     if abs(Y(kk)-Y(kk-1))>.8*pi
%         [ux,inds] = unique(X);
%         uy = Y(inds);
% 
% 
%         % Current value
%         x_me = X(kk);        
%         i_me = find(ux==x_me,1);
% 
%         % Find the previous two phases and compute slope
%         y1 = uy(i_me-1);x1 = ux(i_me-1);
%         y2 = uy(i_me-2);x2 = ux(i_me-2);         
%         dPhi = (y1-y2)/(x1-x2);
% 
%         dXmesign = sign(x_me-x1);
% 
%         Yvec = Y(kk) + [-100:1:150]*2*pi;
% 
%         iLarger = find(Yvec>uy(i_me-1),1);        
%         Yplus = Yvec(iLarger);
%         Yminus = Yvec(iLarger-1);
%         
%         
%         disp(kk);
%         disp(Y(kk));
%         disp(Yplus);
% 
%  
% 
%         if sign(dPhi*dXmesign)==1
%             Y(kk) = Yplus;
%         else
%             Y(kk) = Yminus;
%         end
%     end
end

end
