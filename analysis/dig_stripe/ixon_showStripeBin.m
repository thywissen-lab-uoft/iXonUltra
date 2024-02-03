function [] = ixon_showStripeBin(qgmdata,xVar,opts)

BS = [qgmdata.BinStripe];
P = [qgmdata.Params];
X = [P.(xVar)];

alpha = [BS.ModDepth];
N = [BS.Counts];

goodInds = logical([alpha>.75].*[N>0.5e6]);

P = P(goodInds);
X = X(goodInds);
BS = BS(goodInds);


pos = zeros(length(BS),2);
for kk=1:length(BS)
    s = [BS(kk).Scores];
    c = [BS(kk).Centers];
    
    [s, inds] = sort(s,'descend');
    c = c(inds);
    pos(kk,1) = c(1);
    pos(kk,2) = c(2);
end
hF = figure(9);
clf

subplot(3,2,1);
plot(X,[BS.Phase]/(2*pi),'o');
xlabel(xVar);
if isequal(xVar,'ExecutionDate')
    datetick x
end
ylabel('phase (2\pi)');

subplot(3,2,2);
plot(X,[BS.ModDepth],'o');
xlabel(xVar);
if isequal(xVar,'ExecutionDate')
    datetick x
end
ylabel('ModDepth');
ylim([0 1]);

subplot(3,2,3);
plot(X,pos(:,1),'o');
hold on
% plot(X,pos(:,2),'o');

xlabel(xVar);
if isequal(xVar,'ExecutionDate')
    datetick x
end
ylabel('focus position');

subplot(3,2,4);
plot(X,[BS.Counts],'o');
xlabel(xVar);
if isequal(xVar,'ExecutionDate')
    datetick x
end
ylabel('counts');
% yyaxis right
% plot(X,[BS.RSquareStripe],'o');
subplot(3,2,5);
plot(X,[P.qgm_plane_uwave_frequency_offset],'o');
xlabel(xVar);
if isequal(xVar,'ExecutionDate')
    datetick x
end
ylabel('counts');
end

