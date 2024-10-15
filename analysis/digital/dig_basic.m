function [digdata] = dig_basic(digdata)
    n1 = digdata.n1;
    n2 = digdata.n2;
    [nn1,nn2]=meshgrid(n1,n2);    
     
    % All possible lattice indeces
    N1=nn1(:);
    N2=nn2(:);
    
    
    a1bar = mean(digdata.a1,2);
    a2bar = mean(digdata.a2,2);    
    
    aL = 0.532;
    abar = mean([norm(a1bar) norm(a2bar)]);
    
    % Convert pixel position to site
    px_2_site = @(px_val) px_val/abar;
    
    % Convert pixel position to um position
    px_2_um = @(px_val) aL*(px_val/abar);

            
    for nn=1:length(digdata.FileNames)
        % Get lattice vectors
        a1 = digdata.a1(:,nn);
        a2 = digdata.a2(:,nn);
        A = [a1 a2];                

        % Get lattice phase
        p1 = digdata.p1(:,nn);
        p2 = digdata.p2(:,nn);        
        
        % Get digital data        
        Zdig = digdata.Zdig(:,:,nn);
        
        % Get lattice sites that have atoms
        N1_atom = N1;
        N1_atom(Zdig==0)=[];
        N2_atom = N2;
        N2_atom(Zdig==0)=[];
          
        % Get indeces of every atom
        N=[N1_atom' ; N2_atom'];                      % All points
        p=[p1; p2];
        P = repmat(p,[1 length(N)]);        % Phase vector        
        
        % Compute position of every atom in pixels
        Rn=A*(N+P);         
        
        Xn = Rn(1,:);
        Xc = mean(Xn);
        Xs = std(Xn);
        Yn = Rn(2,:);
        Yc = mean(Yn);
        Ys = std(Yn);        
        
        digdata.RatomSite{nn}=N;
        digdata.Ratom{nn}=Rn;         

        
        digdata.Xc_px(nn) = Xc;
        digdata.Yc_px(nn) = Yc;
        digdata.Xs_px(nn) = Xs;
        digdata.Ys_px(nn) = Ys;   
        
        digdata.Xc_um(nn) = px_2_um(Xc);
        digdata.Yc_um(nn) = px_2_um(Yc);
        digdata.Xs_um(nn) = px_2_um(Xs);
        digdata.Ys_um(nn) = px_2_um(Ys);   
        
        digdata.Xc_site(nn) = px_2_site(Xc);
        digdata.Yc_site(nn) = px_2_site(Yc);
        digdata.Xs_site(nn) = px_2_site(Xs);
        digdata.Ys_site(nn) = px_2_site(Ys);  
                
    end

end

