
function ixondata=ixon_boxCount(ixondata,bgROI)
    fprintf('Performing box count analysis ...');   
    
    if nargin==1
        bgROI=NaN;
    end 
    
    BoxCount=struct;    
    for k=1:size(ixondata.ROI,1)
        ROI=ixondata.ROI(k,:);
        x=ixondata.X(ROI(1):ROI(2));                 % X vector
        y=ixondata.Y(ROI(3):ROI(4));                 % Y vector
        z=double(ixondata.Z(ROI(3):ROI(4),ROI(1):ROI(2)));
        nbg=0;
        
        if nargin==2
            zbg=double(ixondata.Z(bgROI(3):bgROI(4),bgROI(1):bgROI(2)));
            nbg=sum(sum(zbg))/(size(zbg,1)*size(zbg,2)); % count density
        end 
        
        Nraw=sum(sum(z));
        Nbg=nbg*size(z,1)*size(z,2);  
        
        zNoBg=z-nbg;        
        Ncounts=sum(sum(zNoBg));   
        zY=sum(zNoBg,2)';
        zX=sum(zNoBg,1);
               
        % Calculate center of mass
        Xc=sum(zX.*x)/Ncounts;
        Yc=sum(zY.*y)/Ncounts;
        
        % Calculate central second moment/variance and the standard
        % deviation
        X2=sum(zX.*(x-Xc).^2)/Ncounts; % x variance
        Xs=sqrt(X2); % standard deviation X
        Y2=sum(zY.*(y-Yc).^2)/Ncounts; % x variance
        Ys=sqrt(Y2); % standard deviation Y               

        BoxCount(k).Ncounts=Ncounts;    % Number of counts (w/ bkgd removed)
        BoxCount(k).Nraw=Nraw;          % Raw of number of counts
        BoxCount(k).Nbkgd=Nbg;          % Bakcground number of counts
        BoxCount(k).nbkgd=nbg;          % Background counts/px
        BoxCount(k).bgROI=bgROI;        % ROI for calculating bgkd
        BoxCount(k).Xc=Xc;              % X center of mass
        BoxCount(k).Yc=Yc;              % Y center of mass
        BoxCount(k).Xs=Xs;              % X standard deviation
        BoxCount(k).Ys=Ys;              % Y standard deviation
    end    
    ixondata.BoxCount=BoxCount;
    disp('done');
end