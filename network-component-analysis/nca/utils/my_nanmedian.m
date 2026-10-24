function mm = nanmedian(M)

% mm = nanmedian(M)
% medians for each column

mm = nanquantile(M,0.5);
if size(M,1)==2, mm = nanmean(M); end
