function [digdata] = dig_basic(digdata)
%     n1 = digdata.n1; n2 = digdata.n2;       % Site vectors
%     [nn1,nn2]=meshgrid(n1,n2);              % Sites as matrix     
%     N1=nn1(:);N2=nn2(:);                    % lattice indeces    
%     
%     % Lattice vectors
%     a1bar = mean(digdata.a1,2);
%     a2bar = mean(digdata.a2,2);        
%     aL = 0.532; % lattice spacing
%     abar = mean([norm(a1bar) norm(a2bar)]); % average spacing
%     % Convert pixel position to site
%     px_2_site = @(px_val) px_val/abar;    
%     % Convert pixel position to um position
%     px_2_um = @(px_val) aL*(px_val/abar);
    
    a_px = digdata.Lattice_px;
    a_um = digdata.Lattice_um;
            
    for nn=1:length(digdata.FileNames)
        Ratom = digdata.Ratom{nn};
        
        Xc_px   = mean(Ratom(1,:));
        Xc_site = Xc_px/a_px;
        Xc_um   = Xc_site*a_um;
        
        Xs_px   = std(Ratom(1,:));
        Xs_site = Xs_px/a_px;
        Xs_um   = Xs_site*a_um;
        
        Yc_px   = mean(Ratom(2,:));
        Yc_site = Yc_px/a_px;
        Yc_um   = Yc_site*a_um;
        
        Ys_px   = std(Ratom(1,:));
        Ys_site = Ys_px/a_px;
        Ys_um   = Ys_site*a_um;        
          
        digdata.Xc_px(nn) = Xc_px;
        digdata.Xc_site(nn) = Xc_site;
        digdata.Xc_um(nn) = Xc_um;
        digdata.Xs_px(nn) = Xs_px;
        digdata.Xs_site(nn) = Xs_site;
        digdata.Xs_um(nn) = Xs_um;
        
        digdata.Yc_px(nn) = Yc_px;
        digdata.Yc_site(nn) = Yc_site;
        digdata.Yc_um(nn) = Yc_um;
        digdata.Ys_px(nn) = Ys_px;
        digdata.Ys_site(nn) = Ys_site;
        digdata.Ys_um(nn) = Ys_um;

        digdata.nPeakGauss(nn) = 0.5*digdata.Natoms(nn)./(sqrt(2*pi*Xs_site.^2).*sqrt(2*pi*Ys_site.^2));


    
        
%         % Get lattice vectors
%         a1 = digdata.a1(:,nn);
%         a2 = digdata.a2(:,nn);
%         A = [a1 a2];                
% 
%         % Get lattice phase
%         p1 = digdata.p1(:,nn);
%         p2 = digdata.p2(:,nn);        
%         
%         % Get digital data        
%         Zdig = digdata.Zdig(:,:,nn);
%         
%         % Get lattice sites that have atoms
%         N1_atom = N1;
%         N1_atom(Zdig==0)=[];
%         N2_atom = N2;
%         N2_atom(Zdig==0)=[];
%           
%         % Get indeces of every atom
%         N=[N1_atom' ; N2_atom'];                      % All points
%         p=[p1; p2];
%         P = repmat(p,[1 length(N)]);        % Phase vector        
        
%         % Compute position of every atom in pixels
%         Rn=A*(N+P);         
%         
%         Xn = Rn(1,:);
%         Xc = mean(Xn);
%         Xs = std(Xn);
%         Yn = Rn(2,:);
%         Yc = mean(Yn);
%         Ys = std(Yn); 
%         
%         Xc_site = px_2_site(Xc);
%         Yc_site = px_2_site(Yc);
%         Xs_site = px_2_site(Xs);
%         Ys_site = px_2_site(Ys);
        
%         Natoms = digdata.Natoms(nn);      
        
%         npeak = 0.5*Natoms./(sqrt(2*pi*Xs_site.^2).*sqrt(2*pi*Ys_site.^2));
        
%         digdata.RatomSite{nn}=N;
%         digdata.Ratom{nn}=Rn;          
%         digdata.Xc_px(nn) = Xc;
%         digdata.Yc_px(nn) = Yc;
%         digdata.Xs_px(nn) = Xs;
%         digdata.Ys_px(nn) = Ys;           
%         digdata.Xc_um(nn) = px_2_um(Xc);
%         digdata.Yc_um(nn) = px_2_um(Yc);
%         digdata.Xs_um(nn) = px_2_um(Xs);
%         digdata.Ys_um(nn) = px_2_um(Ys);
%         digdata.Xc_site(nn) = Xc_site;
%         digdata.Yc_site(nn) = Yc_site;
%         digdata.Xs_site(nn) = Xs_site;
%         digdata.Ys_site(nn) = Ys_site;          
%         digdata.npeak(nn) = npeak;    
    end

    

end

