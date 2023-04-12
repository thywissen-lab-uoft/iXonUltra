function  [x,fval,exitflag,output]=findLatticePhase(X,Y,Z,k1,k2)


[xx,yy]=meshgrid(X,Y);

zzz = Z(:);
xxx = xx(:);
yyy = yy(:);

myfunc = @(p) ...
    -sum(sin(2*pi*(k1(1)*xxx+k1(2)*k1(2)*yyy)+p(1)).*...
    sin(2*pi*(k2(1)*xxx+k2(2)*k2(2)*yyy)+p(2)).*zzz,'all');
tic
[x,fval,exitflag,output] = fminsearch(myfunc,[0 1]);
toc


    
end

