function [data] = ixon_binnedHistogramFit(data)

if ~isfield(data,'LatticeBin')
    warning('No binned lattice data to analyze')
    return;
end

for kk=1:length(data)
    for nn=1:length(data(kk).LatticeBin)
        z = data(kk).LatticeBin(nn).Zbin;
        % output = bimodalPDFFit(z);
        output = bimodalPDFFit2(z);
        % output = bimodalPDFFit3(z);
        % output = bimodalPDFFit4(z);
        data(kk).LatticeBin(nn).PDFFit = output;
    end
end

end

