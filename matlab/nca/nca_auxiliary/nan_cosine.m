function r = nan_cosine(x,y);

x = repmat(x,size(y,1),1);
x(isnan(y))=nan;
y(isnan(x))=nan;

r = nanmean(x.*y,2) ./ (sqrt(nanmean(x.*x,2)).*sqrt(nanmean(y.*y,2)));
r(~isfinite(r))=-1;

r = 1-r;
