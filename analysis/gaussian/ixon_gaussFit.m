function data=ixon_gaussFit(data,opts)    
    if nargin==1
        opts.doRescale=0;
        opts.doMask=0;
        opts.doRotate=0;
        opts.PCA=[];
        opts.Scale=1;
        opts.Mask=ones(size(data.Z,1),size(data.Z,2));
    end

    GaussFit={};    
    for k=1:size(data.ROI,1)
        ROI=data.ROI(k,:);

        % Acquire the data
        X=double(data.X(ROI(1):ROI(2)));
        Y=double(data.Y(ROI(3):ROI(4)));
        Z=double(data.Z(ROI(3):ROI(4),ROI(1):ROI(2)));       

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
            end
            % Remove data points from fit based on mask
            xx2(~ixon_mask_sc)=[];
            yy2(~ixon_mask_sc)=[];
            Z2(~ixon_mask_sc)=[];
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
            
        else
            % Rotated Gaussian Analysis
            % Very experimental
            % Z is a mxn matrix to which we want to fit an elliptic gaussian            
            
            % Angle of rotation
            theta=atan(opts.PCA.PCA(2,1)/opts.PCA.PCA(1,1));

            % Gaussian Radius
            s1=opts.PCA.Radii(1);
            s2=opts.PCA.Radii(2);
            
            % Center point
            Xc=opts.PCA.Mean(1);
            Yc=opts.PCA.Mean(2);
            
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

            % Perform the fit
            fprintf(' fitting...');
            t1=now;
            [fout,gof,output]=fit([xx2(:) yy2(:)],Z2(:),myfit,opt);
            t2=now;
            disp([' done (' num2str(round((t2-t1)*24*60*60,1)) ' sec.).']);
            GaussFit(k)={fout};
            
        end
    end
    data.GaussFit=GaussFit;

end
