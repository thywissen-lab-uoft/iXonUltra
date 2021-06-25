function ixondata=ixon_fft_boxCount(ixondata,opts)

    fprintf('Performing fft box count analysis ...');    

    if nargin==1
        opts=struct;     
        opts.maskIR=1;
        opts.maskUV=0;
        opts.LMax=100;        
    end
    
    if opts.maskIR
        % Make the mask
        F_IR=1/opts.LMax;
        X=ixondata(1).fft_F;
        Y=ixondata(1).fft_F;
        [xx,yy]=meshgrid(X,Y);
        rr=sqrt(xx.^2+yy.^2);
        IR_mask=rr>F_IR;
    end

    for kk=1:length(ixondata)
        % Initialize structure
        BoxCount=struct;   

        % Get data
        z=ixondata(kk).fft_Z;
        x=ixondata(kk).fft_F;
        y=ixondata(kk).fft_F;
        
        if opts.maskIR
           z=z.*IR_mask; 
        end       

        % Subtract background
        nbg=0;
        Nraw=sum(sum(z));
        Nbg=nbg*size(z,1)*size(z,2);  
        zNoBg=z-nbg;        
        Ncounts=sum(sum(zNoBg));   
        
        % Get sum counts
        zY=sum(zNoBg,2)';
        zX=sum(zNoBg,1);

        zX(zX<0)=0;
        zY(zY<0)=0;

        % Calculate center of mass
        Xc=sum(zX.*x)/Ncounts;
        Yc=sum(zY.*y)/Ncounts;          

        % Calculate central second moment/variance and the standard
        % deviation
        X2=sum(zX.*(x-Xc).^2)/Ncounts; % x variance
        Xs=sqrt(X2); % standard deviation X
        Y2=sum(zY.*(y-Yc).^2)/Ncounts; % x variance
        Ys=sqrt(Y2); % standard deviation Y               

        BoxCount.Ncounts=Ncounts;    % Number of counts (w/ bkgd removed)
        BoxCount.Nraw=Nraw;          % Raw of number of counts
        BoxCount.Nbkgd=Nbg;          % Bakcground number of counts
        BoxCount.nbkgd=nbg;          % Background counts/px
        BoxCount.Xc=Xc;              % X center of mass
        BoxCount.Yc=Yc;              % Y center of mass
        BoxCount.Xs=Xs;              % X standard deviation
        BoxCount.Ys=Ys;              % Y standard deviation
        
        ixondata(kk).fft_BoxCount=BoxCount;
    end 
    
    disp('done');
        
    
end
