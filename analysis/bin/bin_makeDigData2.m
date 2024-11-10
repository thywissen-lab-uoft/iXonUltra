function digdata = bin_makeDigData2(bindata,opts)

if nargin ==1
    opts=struct;
end

if ~isfield(opts,'xVar')
    opts.xVar='ExecutionDate';
end

if ~isfield(opts,'NumSigmaThresh')
    opts.NumSigmaThresh=2.5;
end
    P = [bindata.Params];
    F = [bindata.Flags];
    U = [bindata.Units];
    B = [bindata.LatticeBin];
    PH =[bindata.LatticePhase];
    
    Zdig = zeros(size(bindata(1).LatticeBin(1).Zbin,1),...
        size(bindata(1).LatticeBin(1).Zbin,2),length(bindata));
    n1 = bindata(1).LatticeBin.n1;
    n2 = bindata(1).LatticeBin.n2;
    
    [nn1,nn2]=meshgrid(n1,n2);
    N1_all = nn1(:);
    N2_all = nn2(:);
    
%     Zdig = zeros(numel(n2),numel(n1),length(bindata),length(bindata(1).LatticeBin));
        Zdig = zeros(numel(n2),numel(n1),length(bindata));

    %% Get all Count Centers
    for kk=1:length(bindata)
        for rr=1:length(bindata(kk).LatticeBin)
            count_center(kk,rr) = bindata(kk).LatticeBin(rr).PDF1_Center;
            count_sigma(kk,rr) = bindata(kk).LatticeBin(rr).PDF1_Radius; 
        end
    end 
    
    % Flag outliers
    for nn=1:length(bindata)
        for rr=1:length(bindata(kk).LatticeBin)            
            if (count_center(nn,rr)-median(count_center(:,rr)))>(2*median(count_center(:,rr)))
                count_center(nn,rr)=median(count_center(:,rr));
                count_sigma(nn,rr)=median(count_sigma(:,rr));
            end
        end
    end 
    sigma_norm = count_sigma./count_center;    
    sigma_norm = mean(sigma_norm,1);
    

    Site2PixelsFuncs={};
    Siteatom={};
    Ratom = {};
    Natoms=[];
    lattice_spacing_px=[];
    for nn=1:length(bindata)
        rr=1;
        Zbin = bindata(nn).LatticeBin(rr).Zbin;
        Zscaled = Zbin/count_center(nn,rr);
        Thresh = 1 - opts.NumSigmaThresh*sigma_norm(rr);   
        zdig_this = Zscaled>=Thresh;
        Zdig(:,:,nn) = zdig_this;    
        
        Site2PixelsFuncs{nn}=bindata(nn).LatticeBin(rr).Site2Pixel;
        
        % List of all sites that are occupied
        N1 = N1_all;
        N2 = N2_all;
        N1(~zdig_this(:)) = [];
        N2(~zdig_this(:)) = [];
        
        N1=N1';
        N2= N2';       
        
        Siteatom{nn} =[N1;N2];        
        Ratom{nn}=bindata(nn).LatticeBin(rr).Site2Pixel(N1,N2);
        Natoms(nn) = sum(zdig_this,'all');
        
        a1=norm(bindata(nn).LatticeBin(rr).a1);
        a2=norm(bindata(nn).LatticeBin(rr).a1);
        lattice_spacing_px(nn) = mean([a1 a2]);
    end    
    
    
    
   %% Output

    
    digdata                     = struct;    
    digdata.SourceDirectory     = unique({bindata.SourceDirectory});
    digdata.FileNames           = {bindata.Name}';
    digdata.Threshold           = Thresh;
    digdata.ThresholdType       = 'normalized';
    digdata.xVar                = opts.xVar;
    digdata.X                   = [P.(opts.xVar)];
    digdata.Params              = P;
    digdata.Units               = U;
    digdata.Flags               = F;      
    
    % Actual Digital Data
    digdata.Zdig                = logical(Zdig);
    digdata.n1                  = n1;
    digdata.n2                  = n2;
    digdata.Ratom               = Ratom;
    digdata.Siteatom            = Siteatom;
    digdata.Natoms              = Natoms;
    digdata.Lattice_px          = mean(lattice_spacing_px);
    digdata.Lattice_um          = 1.054/2;

    % Lattice Spacing

    % Basis Information
    digdata.a1                  = [B.a1];       % Lattice Vector 1
    digdata.a2                  = [B.a2];       % Lattice Vector 2    
    digdata.p1                  = [PH.p1];      % Lattice Phase 1
    digdata.p2                  = [PH.p2];      % Lattice Phase 2
    
    
    [digdata] = dig_basic(digdata);