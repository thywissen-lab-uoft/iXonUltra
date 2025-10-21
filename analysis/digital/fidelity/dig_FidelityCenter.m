function digdata = dig_FidelityCenter(digdata)
% Author : CJ Fujiwara
%
% This codes calculates the fidelity of digitized data

% Should there be three figures?  
% (1) Single Shot FIdelity 
% (2) Average image Fidelity
% (3) Summary of Fidelity

if size(digdata.Zdig,4)==1
    warning('cannot compute fidelity on single shot image')
    return;
end

y = digdata.n2;
x = digdata.n1;

R = 30;



[xx,yy]=meshgrid(x,y);
%% Calculate Fidelity
for kk = 1:length(digdata.FileNames)
    z1 = digdata.Zdig(:,:,kk,1);    % Image 1
    z2 = digdata.Zdig(:,:,kk,2);    % Image 2

 
    x1=xx(z1);y1=yy(z1);
    x2=xx(z2);y2=yy(z2);
    
    xc = median(x2);
    yc = median(y2);

    rr = sqrt((xx-xc).^2+(yy-yc).^2);
    rMesh = rr<R;

    z1_sub = z1;z1_sub(~rMesh)=0;
    z2_sub = z2;z2_sub(~rMesh)=0;




    Nlost = sum(z1_sub,'all')-sum(z2_sub,'all');
    Rlost_center = Nlost/sum(z1_sub,'all');


    
    N1 = sum(z1,'all');             % Atom Number 1
    N2 = sum(z2,'all');             % Atom Number 2
    dz = z1 - z2;                   % Differential Image
    Nlost = N1 - N2;                % Number Lost
    Rlost = Nlost/N1;               % Percentage Loss
    Nhop = sum(abs(dz),'all')-Nlost;% Number Hopped
    Rhop = Nhop/N1;                 % Percentage hop

    
    

    digdata.lost_number(kk)       = Nlost;
    digdata.lost_fraction(kk)     = Rlost;
    digdata.lost_fraction_center(kk)  = Rlost_center;
    digdata.hop_number(kk)        = Nhop;
    digdata.hop_fraction(kk)      = Rhop;
end
digdata.FidelityCenterRadius = R;

end
    