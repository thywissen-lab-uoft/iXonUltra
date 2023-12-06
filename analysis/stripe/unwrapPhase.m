
function Y=unwrapPhase(X,Y)
% CF: I think this will be buggy if you have repititions whose noise is
% larger than a phase, but if that happens you're data is too noisy anyway.

% Unwrap the phase with a given X and Y.  Here we shall try to minimize the
% change in slope of the 

[X,inds] = sort(X);
Y = Y(inds);

% Performa a basic unwrap to keep jump less than 2*pi
Y = unwrap(Y);

for kk=3:length(Y)
    if abs(Y(kk)-Y(kk-1))>.8*pi
        [ux,inds] = unique(X);
        uy = Y(inds);


        % Current value
        x_me = X(kk);        
        i_me = find(ux==x_me,1);

        % Find the previous two phases and compute slope
        y1 = uy(i_me-1);x1 = ux(i_me-1);
        y2 = uy(i_me-2);x2 = ux(i_me-2);         
        dPhi = (y1-y2)/(x1-x2);

        dXmesign = sign(x_me-x1);

        Yvec = Y(kk) + [-100:1:150]*2*pi;

        iLarger = find(Yvec>uy(i_me-1),1);        
        Yplus = Yvec(iLarger);
        Yminus = Yvec(iLarger-1);
        
        
        disp(kk);
        disp(Y(kk));
        disp(Yplus);

 

        if sign(dPhi*dXmesign)==1
            Y(kk) = Yplus;
        else
            Y(kk) = Yminus;
        end
    end
end

end
