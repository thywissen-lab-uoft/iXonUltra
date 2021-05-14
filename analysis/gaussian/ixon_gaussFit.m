function ixondata=ixon_gaussFit(ixondata,opts)    
    if nargin==1
        opts.doRescale=0;
        opts.doMask=0;
        opts.doRotate=0;
        opts.Scale=1;
        opts.Mask=ones(size(ixondata.Z,1),size(ixondata.Z,2));
    end

    if ~isfield(ixondata(1),'PCA')
       opts.doRotate=0; 
    end
    
    if opts.doRotate    
        disp(' Performing 2D rotated gaussian fit.');
    else
        disp(' Performing 2D rotated gaussian fit.');
    end

for n=1:length(ixondata)
    GaussFit={};  
    
    fprintf([num2str(n) ' of ' num2str(length(ixondata)) ' :']);

    for k=1:size(ixondata(n).ROI,1)
        ROI=ixondata(n).ROI(k,:);

        % Acquire the data
        X=double(ixondata(n).X(ROI(1):ROI(2)));
        Y=double(ixondata(n).Y(ROI(3):ROI(4)));
        Z=double(ixondata(n).Z(ROI(3):ROI(4),ROI(1):ROI(2)));       

        % Rescale images for fitting speed
        if opts.doRescale
            Z=imresize(Z,opts.Scale);
            X=imresize(X,opts.Scale);
            Y=imresize(Y,opts.Scale);
        end
        
        % Make a mesh grid for fitting
        [xx,yy]=meshgrid(X,Y);    

        % Copy the data
        Z2=Z;xx2=xx;yy2=yy;

        % Remove mask from fit
        if opts.doMask
            % Scale the image mask
            if opts.doRescale                
                ixon_mask_sc=imresize(opts.Mask,opts.Scale);
                xx2(~ixon_mask_sc)=[];
                yy2(~ixon_mask_sc)=[];
                Z2(~ixon_mask_sc)=[];
            else
                ixon_mask=opts.Mask;
                xx2(~ixon_mask)=[];
                yy2(~ixon_mask)=[];
                Z2(~ixon_mask)=[];
            end
            % Remove data points from fit based on mask
            
        end  
        
        if ~opts.doRotate

            % Smooth data, extract peak
            dSmooth=imgaussfilt(Z,0.5);   
            N0=max(max(dSmooth));          

            % Calculate guesses for center and size
            Zx=sum(Z,1);Zy=sum(Z,2)';             % Get X and Y sum profiles
            Nx=sum(Zx);Ny=sum(Zy);                % Get the total number of counts
            Xc=mean(X(Zx>.9*max(X)));           % X center (use >90% SNR)
            Yc=mean(Y(Zy>.9*max(Y)));           % Y center (use >90% SNR)
            Xs=sqrt(sum((X-Xc).^2.*Zx)/Nx); % X standard deviation
            Ys=sqrt(sum((Y-Yc).^2.*Zy)/Ny); % Y standard deviation   
     
            % For ixon hopefully the background is 
            bg=min(min(Z));             
            
            s1=Xs;
            s2=Ys;
            nbg=min(min(Z));
            A=N0;  
             
            gaussme=@(A,Xc,Yc,s1,s2,nbg,xx,yy) A*exp(-( ...
                (1/(2*s1^2))*(xx-Xc).^2 + ...     
                 (1/(2*s2^2))*(yy-Yc).^2))+nbg;        

            myfit=fittype(@(A,Xc,Yc,s1,s2,nbg,xx,yy) gaussme(A,Xc,Yc,s1,s2,nbg,xx,yy),...
                'independent',{'xx','yy'},'coefficients',{'A','Xc','Yc','s1','s2','nbg'});
            opt=fitoptions(myfit);
            opt.StartPoint=[A Xc Yc s1 s2 nbg];
            opt.Upper=[A*1.3 Xc+50 Yc+50 s1*1.5 s2*1.5 1000];
            opt.Lower=[A*0.7 Xc-50 Yc-50 s1*.5 s2*.5 -500];   
            
            % Display initial guess            
            gStr=[' guess (Xc,Yc,s1,s2,A,bg)=(' num2str(round(Xc)) ',' ...
                num2str(round(Yc)) ',' num2str(round(s1)) ',' num2str(round(s2)) ',' ...
                num2str(A,'%.e') ',' num2str(round(nbg)) ')' ];     
            fprintf([gStr ' ... ']);
            
            % Perform the fit
            t1=now;
            [fout,gof,output]=fit([xx2(:) yy2(:)],Z2(:),myfit,opt);
            t2=now;            
            
            % Fit String
            fStr=['(' num2str(round(fout.Xc)) ',' ...
                num2str(round(fout.Yc)) ',' num2str(round(fout.s1)) ',' num2str(round(fout.s2)) ',' ...
                num2str(fout.A,'%.e') ',' num2str(round(fout.nbg)) ')' ];     
            fprintf([' fit ' fStr]);
            
            % dispaly Summary
            disp([' (' num2str(round((t2-t1)*24*60*60,1)) ' sec.).']);
            
            % Append fit to output
            GaussFit(k)={fout}; 
            
        else
            % Rotated Gaussian Analysis
            % Very experimental
            % Z is a mxn matrix to which we want to fit an elliptic gaussian            
            
            PCA=ixondata(n).PCA;
            
            % Angle of rotation     
            theta=atan(PCA.PCA(2,1)/PCA.PCA(1,1));

            % Gaussian Radius
            s1=PCA.Radii(1);
            s2=PCA.Radii(2);
            
            % Center point
            Xc=PCA.Mean(1);
            Yc=PCA.Mean(2);
            % Amplitdue
            A=max(max(imgaussfilt(Z,.5)));
            
            % Background guess
            nbg=min(min(Z));
            
            % https://en.wikipedia.org/wiki/Gaussian_function
            % But we add a minus sign to make it counter clockwise angle
            % When theta=0 s1 is on the x axis
            gaussrot=@(A,Xc,Yc,s1,s2,theta,nbg,xx,yy) A*exp(-( ...
                (cos(theta)^2/(2*s1^2)+sin(theta)^2/(2*s2^2))*(xx-Xc).^2 + ...
                 2*(sin(2*theta)/(4*s1^2) - sin(2*theta)/(4*s2^2))*(xx-Xc).*(yy-Yc) + ...
                 (sin(theta)^2/(2*s1^2)+cos(theta)^2/(2*s2^2))*(yy-Yc).^2))+nbg;           


            myfit=fittype(@(A,Xc,Yc,s1,s2,theta,nbg,xx,yy) gaussrot(A,Xc,Yc,s1,s2,theta,nbg,xx,yy),...
                'independent',{'xx','yy'},'coefficients',{'A','Xc','Yc','s1','s2','theta','nbg'});
            opt=fitoptions(myfit);
            opt.StartPoint=[A Xc Yc s1 s2 theta nbg];
            opt.Upper=[A*1.3 Xc+50 Yc+50 s1*1.5 s2*1.5 theta+5*(pi/180) 1000];
            opt.Lower=[A*0.7 Xc-50 Yc-50 s1*.5 s2*.5 theta-5*(pi/180) -500];

            % Display initial guess  
            gStr=[' guess (Xc,Yc,s1,s2,theta,A,bg)=(' num2str(round(Xc)) ',' ...
                num2str(round(Yc)) ',' num2str(round(s1)) ',' num2str(round(s2)) ',' num2str(round(theta*pi/180,1)) ',' ...
                num2str(A,'%.e') ',' num2str(round(nbg)) ')' ];     
            fprintf([gStr ' ... ']);

            % Perform the fit
            t1=now;
            [fout,gof,output]=fit([xx2(:) yy2(:)],Z2(:),myfit,opt);
            t2=now;            
            
            % Fit String
            fStr=['(' num2str(round(fout.Xc)) ',' ...
                num2str(round(fout.Yc)) ',' num2str(round(fout.s1)) ',' num2str(round(fout.s2)) ',' num2str(round(fout.theta*pi/180,1)) ',' ...
                num2str(fout.A,'%.e') ',' num2str(round(fout.nbg)) ')' ];     
            fprintf([' fit ' fStr]);
            
            % dispaly Summary
            disp([' (' num2str(round((t2-t1)*24*60*60,1)) ' sec.).']);
            
            % Append fit to output
            GaussFit(k)={fout};
            
        end
    end
    ixondata(n).GaussFit=GaussFit;
end

end
