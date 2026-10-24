function mm = nanmean(M,dim)

% mm = nanmean(M)
% means for each column

if ~exist('dim','var'), dim = 1;
  if size(M,1)==1, dim=2; end
end

if dim==2, M=M'; end

for it = 1:size(M,2),
  if sum(isfinite(M(:,it))),
%    M(isfinite(M(:,it)),it)
    mm(it) = mean(M(isfinite(M(:,it)),it));
  else,  mm(it) = nan;
  end
end

if dim==2, mm=mm'; end
