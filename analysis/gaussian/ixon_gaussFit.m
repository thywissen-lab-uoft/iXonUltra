function data=ixon_gaussFit(data,opts)    
    if nargin==1
        opts.doRescale=0;
        opts.doMask=0;
        opts.Scale=1;
        opts.Mask=ones(size(data.Z,1),size(data.Z,2));
    end

    GaussFit={};    
    for k=1:size(data.ROI,1)
        ROI=data.ROI(k,:);

        % Acquire the data
        X=double(data.X(ROI(1):ROI(2)));
        Y=double(data.Y(ROI(3):ROI(4)));
        Zdata=double(data.Z(ROI(3):ROI(4),ROI(1):ROI(2)));       

        % Rescale images for fitting speed
        if opts.doRescale
            Zdata=imresize(Zdata,opts.Scale);
            X=imresize(X,opts.Scale);
            Y=imresize(Y,opts.Scale);
        end

        % Smooth data, extract peak
        dSmooth=imgaussfilt(Zdata,0.5);   
        N0=max(max(dSmooth));          

        % Calculate guesses for center and size
        Zx=sum(Zdata,1);Zy=sum(Zdata,2)';             % Get X and Y sum profiles
        Nx=sum(Zx);Ny=sum(Zy);                % Get the total number of counts
        Xc=mean(X(Zx>.9*max(X)));           % X center (use >90% SNR)
        Yc=mean(Y(Zy>.9*max(Y)));           % Y center (use >90% SNR)
        Xs=sqrt(sum((X-Xc).^2.*Zx)/Nx); % X standard deviation
        Ys=sqrt(sum((Y-Yc).^2.*Zy)/Ny); % Y standard deviation       
  
        % Make a mesh grid for fitting
        [xx,yy]=meshgrid(X,Y);

        % Copy the data
        Z2=Zdata;xx2=xx;yy2=yy;
        
        % Remove mask from fit
        if opts.doMask
            % Scale the image mask
            if opts.doRescale                
                ixon_mask_sc=imresize(opts.Mask,opts.Scale);
            end
            % Remove data points from fit based on mask
            xx2(~ixon_mask_sc)=[];
            yy2(~ixon_mask_sc)=[];
            Z2(~ixon_mask_sc)=[];
        end      

        % For ixon hopefully the background is 
        bg=min(min(Zdata)); 

        % Create fit object
        myfit=fittype('A*exp(-(xx-Xc).^2./(2*Xs^2)).*exp(-(yy-Yc).^2./(2*Ys^2))+nbg',...
            'independent',{'xx','yy'},'coefficients',{'A','Xc','Xs','Yc','Ys','nbg'});
        opt=fitoptions(myfit);
        opt.StartPoint=[N0 Xc Xs Yc Ys bg];
        opt.Lower=[N0/10 10 1 10 1 -1];
        opt.Upper=[5*N0 max(X) range(X) max(Y) range(Y) N0];
        opt.Weights=[];
        

        % Check that the upper and lower bounds make sense
        badInds=opt.Upper<opt.Lower;
        if sum(badInds)
            warning(['Generated lower bounds for gaussian fit exceed the upper ' ...
                'bounds for ' num2str(sum(badInds)) ' parameters. This ' ...
                'may be caused by no atoms.']);
            opt.Lower=[0 0 0 0 0 0];
            opt.Upper=[];
            opt.StartPoint=[100 mean(Dx) range(Dx)/10 mean(Dy) range(Dy)/10 ...
                10];
        end

        % Display initial guess
        str1=['(Xc0,Yc0)=(' num2str(round(Xc)) ',' num2str(round(Yc)) ');'];
        str2=['(Xs0,Ys0)=(' num2str(round(Xs)) ',' num2str(round(Ys)) ')'];
        fprintf([str1 str2 ';']);

        % Perform the fit
        fprintf(' fitting...');
        t1=now;
        [fout,gof,output]=fit([xx2(:) yy2(:)],Z2(:),myfit,opt);
        t2=now;
        disp([' done (' num2str(round((t2-t1)*24*60*60,1)) ' sec.).']);
        GaussFit(k)={fout};
    end
    data.GaussFit=GaussFit;

end