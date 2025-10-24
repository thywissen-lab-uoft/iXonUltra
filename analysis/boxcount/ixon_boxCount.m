function ixondata=ixon_boxCount(ixondata,bgROI)

    fprintf('Performing box count analysis ...');    
    if nargin==1
        disp(' No background ROI provided, will assume background of zero.');
        bgROI=[];          
    else
        disp([' Using background counts from ROI = [' ...
            num2str(bgROI) ']']);        
    end    
    
    for kk=1:length(ixondata)

        BoxCount=struct;    
        for ll=1:size(ixondata(kk).Z,3)
            for k=1:size(ixondata(kk).ROI,1)
                ROI=ixondata(kk).ROI(k,:);
              
                ix_1 = find(ixondata(kk).X>=ROI(1),1);
                ix_2 = find(ixondata(kk).X>=ROI(2),1);
                iy_1 = find(ixondata(kk).Y>=ROI(3),1);
                iy_2 = find(ixondata(kk).Y>=ROI(4),1);            
                x = ixondata(kk).X(ix_1:ix_2);
                y = ixondata(kk).Y(iy_1:iy_2);    
                z = ixondata(kk).Z(iy_1:iy_2,ix_1:ix_2,ll);  
                nbg=0;  
                Npeak = max(z,[],'all');            
                Nraw=sum(sum(z));
                Nbg=nbg*size(z,1)*size(z,2);  
                zNoBg=z-nbg;        
                Ncounts=sum(sum(zNoBg));   
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
    
                BoxCount(ll,k).Ncounts=Ncounts;    % Number of counts (w/ bkgd removed)
                BoxCount(ll,k).Npeak=Npeak;        % Peak Counts (no bkgd removed)
                BoxCount(ll,k).Nraw=Nraw;          % Raw of number of counts
                BoxCount(ll,k).Nbkgd=Nbg;          % Bakcground number of counts
                BoxCount(ll,k).nbkgd=nbg;          % Background counts/px
                BoxCount(ll,k).bgROI=bgROI;        % ROI for calculating bgkd
                BoxCount(ll,k).Xc=Xc;              % X center of mass
                BoxCount(ll,k).Yc=Yc;              % Y center of mass
                BoxCount(ll,k).Xs=Xs;              % X standard deviation
                BoxCount(ll,k).Ys=Ys;              % Y standard deviation
            end  
        end
        ixondata(kk).BoxCount=BoxCount;
    end
end
