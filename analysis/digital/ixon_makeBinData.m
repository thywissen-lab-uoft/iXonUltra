function bindata = ixon_makeBinData(ixondata)
% Initialize outputs
bindata = struct;

for nn = 1:length(ixondata)
    bindata(nn).SourceDirectory = fileparts(saveOpts.saveDir);
    bindata(nn).Date = ixondata(nn).Date;
    bindata(nn).Name = ixondata(nn).Name;
    bindata(nn).Params = ixondata(nn).Params;
    bindata(nn).Units = ixondata(nn).Units;
    bindata(nn).Flags = ixondata(nn).Flags;
    bindata(nn).CameraInformation = ixondata(nn).CameraInformation;
    bindata(nn).AcquisitionInformation = ixondata(nn).AcquisitionInformation;
    bindata(nn).AcquisitionDescription = ixondata(nn).AcquisitionDescription;
    bindata(nn).ROI = ixondata(nn).ROI;
    bindata(nn).ProcessOptions = ixondata(nn).ProcessOptions;  
    bindata(nn).LatticeK = ixondata(nn).LatticeK;
    bindata(nn).LatticePhase = ixondata(nn).LatticePhase;
    bindata(nn).LatticeBin = ixondata(nn).LatticeBin;
end

end

